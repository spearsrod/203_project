from gym_test import gym
import pyglet
from pyglet import gl
import numpy as np
import time

FPS = 50


def get_path_curvature(xy_line):
	k = np.zeros(xy_line.shape[0])
	for idx in range(1,xy_line.shape[0] - 1):
		x_cur = xy_line[idx,0]
		y_cur = xy_line[idx,1]
		x_prev = xy_line[idx-1,0]
		y_prev = xy_line[idx-1,1]
		x_next = xy_line[idx+1,0]
		y_next = xy_line[idx+1,1]

		A = np.array([[x_cur, y_cur, 1],[x_prev, y_prev, 1],[x_next, y_next, 1]])
		det = np.linalg.det(A)
		if(det == 0):
			k[idx] =0
		else:
			b = np.array([[-(x_cur**2) - (y_cur**2)],[-(x_prev**2) - (y_prev**2)],[-(x_next**2) - (y_next**2)]])

			coefs = np.linalg.solve(A, b)

			R = np.sqrt(np.square(coefs[0])/4 + np.square(coefs[1])/4 - coefs[2])
			k[idx] = 1.0/R
	return k

def extract_center_line(env):
	track_array = np.asarray(env.track)
	track_center = track_array[:,2:]
	return track_center

def interpolate_path(path, N):
	x_path = path[:,0]
	y_path = path[:,1]

	t = np.arange(0,path.shape[0])
	t_upsampled = np.linspace(0, path.shape[0], N)

	x_interp = np.interp(t_upsampled, t, x_path)
	y_interp = np.interp(t_upsampled, t, y_path)

	interp_path = np.zeros((N,2))
	interp_path[:,0] = x_interp
	interp_path[:,1] = y_interp
	return interp_path

def cartesian_to_path(xy_line):
	s = np.zeros(xy_line.shape[0])
	for idx in range(1,xy_line.shape[0]):
		x_cur = xy_line[idx,0]
		y_cur = xy_line[idx,1]
		x_prev = xy_line[idx-1,0]
		y_prev = xy_line[idx-1,1]
		dx = x_cur - x_prev
		dy = y_cur - y_prev
		s[idx] = np.sqrt(dx**2 + dy**2) + s[idx - 1]

	k = get_path_curvature(xy_line)
	return s, k

def get_car_position(env):
	return np.asarray(env.car.hull.position)

def render_center_line(env):
	# track_center = extract_center_line(env)
	# gl.glBegin(gl.GL_QUADS)
	# gl.glColor4f(1.0,0,0,1.0)
	# for idx in range(0,track_center.shape[0]):
	# 	cur_point = track_center[idx,:]
	# 	size = 10
	# 	gl.glVertex3f(cur_point[0] + size,cur_point[1]+ size,0)
	# 	gl.glVertex3f(cur_point[0] + size,cur_point[1]- size,0)
	# 	gl.glVertex3f(cur_point[0] - size,cur_point[1]+ size,0)
	# 	gl.glVertex3f(cur_point[0] - size,cur_point[1]- size,0)
	# gl.glEnd()

	pos = get_car_position(env)
	gl.glColor4f(1,0,0,1)
	gl.glBegin(gl.GL_POINTS)
	gl.glVertex3f(pos[0] + 4, pos[1] + 2, 0)
	gl.glEnd()

def find_nearest_path(pos_cur, path):
	distance = np.sqrt(np.square(pos_cur[0] - path[:,0]) + np.square(pos_cur[1] - path[:,1]))
	idx = np.argmin(distance)
	sign = np.sign(pos_cur[1] - path[idx,0])

	return idx, sign*distance[idx]

def update_dynamics(env, path_s, path_k, path_state):
	pos_cur = np.array([np.asarray(env.car.hull.position)[0],np.asarray(env.car.hull.position)[1],env.car.hull.angle])
	dt = 1.0/FPS
	Ux = env.car.hull.vehicle_velocity[0]
	Uy = env.car.hull.vehicle_velocity[1]

	s = path_state[0]
	e = path_state[1]
	d_psi = path_state[2]

	path_idx = (np.abs(s - path_s)).argmin()
	K = path_k[path_idx]
	#r = (pos_cur[2] - pos_prev[2])/dt
	r = env.car.hull.yaw_rate

	

	s_dot = (1/(1-e*K))*(Ux*np.cos(d_psi) - Uy*np.sin(d_psi))
	e_dot = -Uy*np.cos(d_psi) - Ux*np.sin(d_psi)
	d_psi_dot = r - K*s_dot

	s = s + s_dot*dt
	e = e + e_dot*dt
	d_psi = d_psi + d_psi_dot*dt

	return np.array([s, e, d_psi]), pos_cur

def lookahead_controller(env, path, path_state, car_state, speed_profile, prev_max):
	car = env.car
	Kla = 3500
	xla = 15

	e = path_state[1]
	d_psi = path_state[2]
	path_idx, e = find_nearest_path(car_state, path)

	Ux_des = speed_profile[path_idx]

	Kd = 0.00000001
	Ux = car.hull.vehicle_velocity[0]
	lat_force = 0
	Fxtotal = Kd*(Ux_des - Ux)
	print(Fxtotal)
	if(Fxtotal > prev_max):
		lat_force = 1
		new_max = Fxtotal
	elif(Fxtotal < 0):
		lat_force = -1
		new_max = prev_max
	else:
		print('does this happen?')
		lat_force = Fxtotal/prev_max/4.0
		print(lat_force)
		new_max = prev_max


	K = path[path_idx,1]
	dif = env.car.hull.angle - np.asarray(env.track)[0,1]
	dif_sign = np.sign(np.mod(dif,np.pi/2.0) - np.mod(dif,np.pi))
	if(dif_sign ==0):
		dif_sign = 1

	# dif_sign*np.mod(dif,np.pi/2.0)


	Caf = 275000
	delta = -Kla*(e + xla*d_psi)/Caf
	return -delta, lat_force, new_max

def generate_speed_profile(path_s, path_k, env):
	F_max = env.car.friction*env.car.m*9.81
	a_max = (F_max/env.car.m)

	profile = np.zeros(path_k.size)

	lookahead = 2
	for idx in range(0,path_k.size):
		if(idx + lookahead < path_k.size):
			average_k = np.average(path_k[idx:idx+lookahead])
		else:
			average_k = np.average(path_k[idx:])

		max_Ux = np.sqrt(1.0/(average_k + 0.0001)*a_max)
		profile[idx] = np.minimum(max_Ux
	return profile



def main():
	env = gym.make('CarRacing-v0')

	
	env.reset()

	center_line = extract_center_line(env)

	interp_center_line = interpolate_path(center_line, center_line.shape[0]*2)

	interp_s, interp_k = cartesian_to_path(interp_center_line)
	path_state = np.array([interp_s[0], 0, 0])
	speed_profile = generate_speed_profile(interp_s, interp_k, env)
	new_max = 0
	for idx in range(1000):
		env.render()

		path_state, car_state = update_dynamics(env, interp_s, interp_k, path_state)
		delta, Fx, new_max = lookahead_controller(env, interp_center_line, path_state, car_state, speed_profile, new_max)
		if(Fx > 0):
			action = [delta, 1, 0]
		else:
			action = [delta, 0, 1]
		# action = [delta, 0.8, 0]
		# if(idx > 150):
		# 	action = [delta, 0, 0.6]
		# if(idx > 300):
		# 	action = [delta, 0.4, 0]
		# if(idx > 400):
		# 	action = [delta, 0, 0]
		# print(action[1])
		# print('hi')
		# print(env.car.hull.position)
		# print(env.track[0][1:4])
		env.step(action)
	time.sleep(20)
	env.close()


if __name__ == '__main__':
	main()