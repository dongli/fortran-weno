#!/usr/bin/env python3

from sympy import *
init_printing()

x, y, dx, dy = symbols('x y \Delta{x} \Delta{y}')

def reconstruct_polynomial_1d(m):
	'''
		Reconstruct 1D polynomial from averaged values on the stencil.

		Input:
			m - Stencil size (e.g. 5 for WENO5)

		Output:
			f(x) - Reconstruced polynomial
	'''
	a = Matrix(symbols(f'a(0:{m})'))
	x, xi, dx = symbols('x x_i \Delta{x}')
	poly_basis = Matrix([x**i for i in range(m)])

	f = sum(HadamardProduct(a, poly_basis, evaluate=True))

	F = simplify(integrate(f, (x, (xi - dx / 2, xi + dx / 2))) / dx)
	for k in range(m): F = F.collect(a[k])
	F = F.subs(xi, x)

	tmp = []
	for i in range(m):
		tmp.append(symbols('F_{' + str(i) + '}'))
	Fi = Matrix(tmp)

	tmp = []
	for cx in range(-int((m-1)/2), int((m-1)/2)+1):
		g = expand(simplify(F.subs(x, cx * dx)))
		tmp.append([])
		for k in range(m):
			tmp[-1].append(g.coeff(a[k]))
	A = Matrix(tmp)

	ak = list(simplify(linsolve(A * a - Fi, list(a))))[0]

	for k in range(m): f = f.subs(a[k], ak[k])
	f = expand(f)

	return f

def reconstruct_polynomial_2d(m, n, Fi=None, mask=None):
	'''
		Reconstruct 2D polynomial from averaged values on the stencils.

		Input:
			m - Stencil size along x-axis
			n - Stencil size along y-axis

		Output:
			f(x) - Reconstructed polynomial
	'''
	tmp = []
	for j in range(m):
		for i in range(n):
			if mask == None or mask[j,i] == 0:
				tmp.append(symbols('a_{' + str(i) + str(j) + '}'))
	a = Matrix(tmp)

	x, y, xi, yi, dx, dy = symbols('x y x_i y_i \Delta{x} \Delta{y}')
	tmp = []
	for j in range(m):
		for i in range(n):
			if mask == None or mask[j,i] == 0:
				tmp.append(x**i * y**j)
	p = Matrix(tmp)

	f = sum(HadamardProduct(a, p, evaluate=True))

	F = simplify(integrate(integrate(f, (x, (xi - dx / 2, xi + dx / 2))), (y, (yi - dy / 2, yi + dy / 2))) / (dx * dy))
	for k in range(len(a)):
		F = F.collect(a[k])
	F = expand(F.subs(xi, x).subs(yi, y))

	tmp = []
	j = 0
	for cy in range(-int((n-1)/2), int((n-1)/2)+1):
		i = 0
		for cx in range(-int((m-1)/2), int((m-1)/2)+1):
			if mask == None or mask[j,i] == 0:
				g = expand(simplify(F.subs(x, cx * dx).subs(y, cy * dy)))
				tmp.append([])
				for k in range(len(a)):
					tmp[-1].append(g.coeff(a[k]))
			i += 1
		j += 1
	A = Matrix(tmp)

	if not Fi:
		tmp = []
		for j in range(m):
			for i in range(n):
				if mask == None or mask[j,i] == 0:
					tmp.append(symbols('F_{' + str(i) + str(j) + '}'))
		Fi = Matrix(tmp)
	elif mask != None:
		tmp = []
		k = 0
		for j in range(m):
			for i in range(n):
				if mask == None or mask[j,i] == 0:
					tmp.append(Fi[k])
				k += 1
		Fi = Matrix(tmp)

	a = Matrix(list(a))
	ak = list(simplify(linsolve(A * a - Fi, list(a))))[0]

	for k in range(len(a)): f = f.subs(a[k], ak[k])
	f = expand(f)

	return f

f_5 = reconstruct_polynomial_1d(5)

print(list(f_5.subs(x, dx / 2).as_coefficients_dict().values()))

f_3x3 = reconstruct_polynomial_2d(3, 3)

print(list(f_3x3.subs(x, dx / 2).subs(y, dy / 2).as_coefficients_dict().values()))
