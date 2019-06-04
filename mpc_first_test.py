import cvxpy as cp
import numpy as np
from gym_test import gym


SCALE = 6.0
TRACK_WIDTH = 40/SCALE
FPS         = 50
Iz = 2235
m = 1648
a = 1.234
b = 1.234
Ca_f = 275000
Ca_r = 265000


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

def direct_optimization(N, dt, path_state, Ux, path_idx, ds):
	print('check')
	print(path_state)
	s0 = path_state[0]
	s_dot0 = path_state[1]
	e0 = path_state[2]
	e_dot0 = path_state[3]
	dpsi0 = path_state[4]
	dpsi_dot0 = path_state[5]

	est_idx = Ux*N*dt/ds
	split_idx = est_idx/N
	K = np.rint(split_idx*np.arange(0,N))
	K_dot = np.zeros(N)
	for idx in range(1,N):
		K_dot[idx] = K[idx] - K[idx - 1]


	delta = cp.Variable((N,1))
	s = cp.Variable((N,1))
	s_dot = cp.Variable((N,1))
	e = cp.Variable((N,1))
	e_dot = cp.Variable((N,1))
	dpsi = cp.Variable((N,1))
	dpsi_dot = cp.Variable((N,1))

	f0 = -s[-1]
	objective = cp.Minimize(f0)
	#subject to x_i+1 = x_i + hAx_i + hBu_i + hC
	constraints = [s[0] == s0, e[0] == e0, dpsi[0]==dpsi0, dpsi_dot[0] == dpsi_dot0, e_dot[0] == e_dot0, s_dot[0] == Ux]
	constraints += [delta[0] <= np.radians(25), delta[0] >= np.radians(-25)]
	for idx in range(1, N):
		constraints += [s[idx] == s[idx - 1] + s_dot[idx - 1]*dt]
		constraints += [s_dot[idx] == Ux]
		constraints += [e[idx] == e[idx - 1] + e_dot[idx - 1]*dt]
		constraints += [e_dot[idx] == e_dot[idx - 1] + dt*(-(Ca_f + Ca_r)/(m*Ux)*e_dot[idx - 1] + (Ca_f + Ca_r)/m*dpsi[idx - 1] + (b*Ca_r - a*Ca_f)/(m*Ux)*dpsi_dot[idx - 1]) + dt*(Ca_f/m)*delta[idx-1] + dt*(-K[idx - 1]*(Ux**2 - (b*Ca_r - a*Ca_f)/m))]
		constraints += [dpsi[idx] == dpsi[idx - 1] + dpsi_dot[idx - 1]*dt]
		constraints += [dpsi_dot[idx] == dpsi_dot[idx - 1] + dt*((b*Ca_r - a*Ca_f)/(Iz*Ux)*e_dot[idx - 1] + (a*Ca_f - b*Ca_r)/Iz*dpsi[idx - 1] - (a**2*Ca_f + b**2*Ca_r)/(Iz*Ux)*dpsi_dot[idx - 1]) + dt*a*Ca_f/Iz*delta[idx - 1] + dt*(-K[idx - 1]*(a**2*Ca_f + b**2*Ca_r)/Iz - K_dot[idx - 1]*Ux)]
		constraints += [delta[idx] <= np.radians(25), delta[idx] >= np.radians(-25)]
		constraints += [e[idx] <= TRACK_WIDTH/2.0, e[idx] >= -TRACK_WIDTH/2.0]


	prob = cp.Problem(objective, constraints)
	result = prob.solve()

	return -delta.value[0][0]

def find_nearest_path(pos_cur, path):
	distance = np.sqrt(np.square(pos_cur[0] - path[:,0]) + np.square(pos_cur[1] - path[:,1]))
	idx = np.argmin(distance)
	sign = np.sign(pos_cur[1] - path[idx,0])

	return idx, sign*distance[idx]

def update_path(env, center_line, prev_state, s, k):
	pos_cur = np.array([env.car.hull.position[0], env.car.hull.position[1]])
	idx, e = find_nearest_path(pos_cur, center_line)

	next_state = np.zeros(6)
	next_state[0] = s[idx]
	next_state[1] = s[idx] - prev_state[0]
	next_state[2] = e
	next_state[3] = e - prev_state[2]
	psi_path = np.arctan2(center_line[idx+1,1] - center_line[idx,1], center_line[idx+1,0] - center_line[idx,0])
	next_state[4] = env.car.hull.angle - psi_path
	next_state[5] = next_state[4] - prev_state[4]
	return next_state, idx



def main():
	env = gym.make('CarRacing-v0')

	
	env.reset()

	N = 10
	dt = 1.0/FPS
	env.render()


	center_line = extract_center_line(env)

	interp_center_line = interpolate_path(center_line, center_line.shape[0]*1)

	interp_s, interp_k = cartesian_to_path(interp_center_line)

	Ux = 10.0
	path_state = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
	path_idx = 0
	ds = interp_s[1] - interp_s[0]

	#delta = direct_optimization(N, dt, path_state, Ux, path_idx, ds)
	print(env.car.hull.position)
	print(ds)
	for idx in range(1000):
		env.render()

		delta = direct_optimization(N, dt, path_state, Ux, path_idx, ds)
		ENGINE_POWER = env.car.m*1.0*9.81
		action = [0.0, 1.0*ENGINE_POWER, 0.0]
		env.step(action)
		print(env.car.hull.position)
		path_state, path_idx = update_path(env, interp_center_line, path_state, interp_s, interp_k)
	time.sleep(20)
	env.close()


if __name__ == '__main__':
	main()