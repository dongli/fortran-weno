#!/usr/bin/env python3

from sympy import *
init_printing()

x, y, dx, dy = symbols('x y \Delta{x} \Delta{y}')

def reconstruct_polynomial_1d(m, Fi=None, mask=None):
	'''
		Reconstruct 1D polynomial from averaged values on the stencil.

		Input:
			m    - Stencil size (e.g. 5 for WENO5)
			Fi   - Symbols for averaged values
			mask - Mask missing cells

		Output:
			f(x) - Reconstruced polynomial
	'''
	x, xi, dx = symbols('x x_i \Delta{x}')

	# Set symbols for polynomial coefficients.
	tmp = []
	for i in range(m):
		if mask == None or mask[i] == 0: tmp.append(symbols(f'a_{i}'))
	a = Matrix(tmp)

	# Set polynomial basis which are monomials.
	tmp = []
	if mask != None and mask[0] == 1:
		for i in range(m):
			if mask[m-i-1] == 0: tmp.append(x**i)
	else:
		for i in range(m):
			if mask == None or mask[i] == 0: tmp.append(x**i)
	p = Matrix(tmp)

	# Set polynomial function.
	f = sum(HadamardProduct(a, p, evaluate=True))

	# Integrate polynomial function within a cell.
	F = simplify(integrate(f, (x, (xi - dx / 2, xi + dx / 2))) / dx)
	for k in range(len(a)): F = F.collect(a[k])
	F = F.subs(xi, x)

	# Calculate coefficient matrix of integrated function on each cell centers.
	tmp = []; i = 0
	if m % 2 == 0:
		cell_coefs = range(-int(m/2), int(m/2))
	else:
		cell_coefs = range(-int((m-1)/2), int((m-1)/2)+1)
	for cx in cell_coefs:
		if mask == None or mask[i] == 0:
			g = expand(simplify(F.subs(x, cx * dx)))
			tmp.append([])
			for k in range(len(a)): tmp[-1].append(g.coeff(a[k]))
		i += 1
	A = Matrix(tmp)

	tmp = []
	if Fi == None:
		for i in range(m):
			if mask == None or mask[i] == 0: tmp.append(symbols(f'F_{i}'))
	else:
		for i in range(m):
			if mask == None or mask[i] == 0: tmp.append(Fi[i])
	Fi = Matrix(tmp)

	# Solve the equations under the conditions that f should recover the cell averaged values.
	ak = list(simplify(linsolve(A * a - Fi, list(a))))[0]

	for k in range(len(a)): f = f.subs(a[k], ak[k])
	f = expand(f)

	return f

def calc_ideal_weights_1d(m, x0, mask):
	dx = symbols('\Delta{x}')

	# The second element is the relative x coordinate in each sub-stencils.
	#  _________________________________________
	# |     |     |     |     |     |     |     |
	# |     |     |     |  0 x0     |     |     |
	# |_____|_____|_____|_____|_____|_____|_____|
	#
	#  _______________________
	# |     |     |     |     |
	# |     |     |  0  |   x0+dx
	# |_____|_____|_____|_____|
	#
	#        _______________________
	#       |     |     |     |     |
	#       |     |     |  0 x0     |
	#       |_____|_____|_____|_____|
	#
	#              _______________________
	#             |     |     |     |     |
	#             |     |  x0-dx 0  |     |
	#             |_____|_____|_____|_____|
	#
	#                    _______________________
	#                   |     |     |     |     |
	#                   |  x0-2dx   |  0  |     |
	#                   |_____|_____|_____|_____|
	#
	subs = {
		5: (3, (x0 + dx, x0, x0 - dx)),
		7: (4, (x0 + dx, x0, x0 - dx, x0 - 2 * dx))
	}

	# Set symbols for cell averaged values.
	tmp = []
	for i in range(m): tmp.append(symbols(f'F_{i}'))
	Fi = Matrix(tmp)
	list_Fi = list(Fi)

	# Calculate weights for sub-stencils.
	A = zeros(m, subs[m][0])
	i = 0
	for local_x in subs[m][1]:
		local_Fi = Fi[i:i+subs[m][0]]
		local_mask = mask[i:i+subs[m][0]]
		f = reconstruct_polynomial_1d(subs[m][0], Fi=Matrix(local_Fi), mask=local_mask)
		tmp = f.subs(x, local_x).as_coefficients_dict()
		for key in local_Fi: A[list_Fi.index(key),i] = tmp[key]
		i += 1

	# Calculate weights for full stencil.
	f = reconstruct_polynomial_1d(m, mask=mask)
	tmp = f.subs(x, x0).as_coefficients_dict()
	b = zeros(m, 1)
	for k in range(len(b)): b[k] = tmp[list_Fi[k]]

	# Calculate ideal weights.
	w = (transpose(A) * A)**-1 * transpose(A) * b

	# Check result.
	if sum(A * w - b) != 0:
		print('Failed to calculate ideal weights')
		return None

	return w

def reconstruct_polynomial_2d(m, n, Fi=None, mask=None):
	'''
		Reconstruct 2D polynomial from averaged values on the stencil.

		Input:
			m    - Stencil size along x-axis
			n    - Stencil size along y-axis
			Fi   - Symbols for averaged values
			mask - Mask missing cells

		Output:
			f(x) - Reconstructed polynomial
	'''
	x, y, xi, yi, dx, dy = symbols('x y x_i y_i \Delta{x} \Delta{y}')

	# Set symbols for polynomial coefficients.
	tmp = []
	for j in range(n):
		for i in range(m):
			if mask == None or mask[j,i] == 0: tmp.append(symbols('a_{' + str(i) + str(j) + '}'))
	a = Matrix(tmp)

	# The mask used to filter monomials cab be different, the ones should always at the right-top corners.
	tmp = []
	if mask != None and mask[0,0] == 1:
		for j in range(n):
			for i in range(m):
				if mask[n-j-1,m-i-1] == 0: tmp.append(x**i * y**j)
	elif mask != None and mask[-1,0] == 1:
		for j in range(n):
			for i in range(m):
				if mask[j,m-i-1] == 0: tmp.append(x**i * y**j)
	elif mask != None and mask[0,-1] == 1:
		for j in range(n):
			for i in range(m):
				if mask[n-j-1,i] == 0: tmp.append(x**i * y**j)
	else:
		for j in range(n):
			for i in range(m):
				if mask == None or mask[j,i] == 0: tmp.append(x**i * y**j)
	p = Matrix(tmp)

	# Set polynomial function.
	f = sum(HadamardProduct(a, p, evaluate=True))

	# Integrate polynomial function within a cell.
	F = simplify(integrate(integrate(f, (x, (xi - dx / 2, xi + dx / 2))), (y, (yi - dy / 2, yi + dy / 2))) / (dx * dy))
	for k in range(len(a)): F = F.collect(a[k])
	F = expand(F.subs(xi, x).subs(yi, y))

	# Calculate coefficient matrix of integrated function on each cell centers.
	tmp = []; j = 0
	for cy in range(-int((n-1)/2), int((n-1)/2)+1):
		i = 0
		for cx in range(-int((m-1)/2), int((m-1)/2)+1):
			if mask == None or mask[j,i] == 0:
				g = expand(simplify(F.subs(x, cx * dx).subs(y, cy * dy)))
				tmp.append([])
				for k in range(len(a)): tmp[-1].append(g.coeff(a[k]))
			i += 1
		j += 1
	A = Matrix(tmp)

	tmp = []
	if not Fi:
		for j in range(n):
			for i in range(m):
				if mask == None or mask[j,i] == 0:
					tmp.append(symbols('F_{' + str(i) + str(j) + '}'))
	else:
		k = 0
		for j in range(n):
			for i in range(m):
				if mask == None or mask[j,i] == 0:
					tmp.append(Fi[k])
				k += 1
	Fi = Matrix(tmp)

	# Solve the equations under the conditions that f should recover the cell averaged values.
	ak = list(simplify(linsolve(A * a - Fi, list(a))))[0]

	for k in range(len(a)): f = f.subs(a[k], ak[k])
	f = expand(f)

	return f

def calc_ideal_weights_2d(m, x0, y0, mask):
	dx, dy = symbols('\Delta{x} \Delta{y}')

	subs = {
		5: (3, 3, (x0 + dx, x0, x0 - dx), (y0 + dy, y0, y0 - dy)),
	}

	tmp = []
	for j in range(m):
		tmp.append([])
		for i in range(m):
			tmp[-1].append(symbols('F_{' + str(i) + str(j) + '}'))
	Fi = Matrix(tmp)
	list_Fi = list(Fi)

	A = zeros(m*m, subs[m][0]*subs[m][1])
	k = 0; j = 0
	for local_y in subs[m][3]:
		i = 0
		for local_x in subs[m][2]:
			local_Fi = flatten(Fi[j:j+subs[m][1],i:i+subs[m][0]])
			local_mask = mask[j:j+subs[m][1],i:i+subs[m][0]]
			f = reconstruct_polynomial_2d(subs[m][0], subs[m][1], Fi=Matrix(local_Fi), mask=local_mask)
			tmp = f.subs(x, local_x).subs(y, local_y).as_coefficients_dict()
			for key in local_Fi:
				A[list_Fi.index(key),k] = tmp[key]
			i += 1
			k += 1
		j += 1

	f = reconstruct_polynomial_2d(m, m, mask=mask)
	tmp = f.subs(x, x0).subs(y, y0).as_coefficients_dict()
	b = zeros(m*m, 1)
	for k in range(len(b)):
		b[k] = tmp[list_Fi[k]]

	w = (transpose(A) * A)**-1 * transpose(A) * b

	if sum(A * w - b) != 0:
		print('Failed to calculate ideal weights')
		return None

	return w

f_5 = reconstruct_polynomial_1d(5)

print(list(f_5.subs(x, dx / 2).as_coefficients_dict().values()))

f_3x3 = reconstruct_polynomial_2d(3, 3)

print(list(f_3x3.subs(x, dx / 2).subs(y, dy / 2).as_coefficients_dict().values()))
