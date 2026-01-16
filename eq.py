from sympy import symbols, Matrix, Eq, solve

gamma, omega, ux, uy, uz = symbols('g o ux uy uz')
tau, x, y, z = symbols('tau x y z')
tau_pt, x_pt, y_pt, z_pt = symbols(" tau' x' y' z' ")

sigma_x = Matrix([[0,1],
                  [1,0]])
sigma_y = Matrix([[0,-1j],
                  [1j,0]])
sigma_z = Matrix([[1,0],
                  [0,-1]])
sigma_plus = (sigma_x + 1j*sigma_y)/2
sigma_moins = (sigma_x - 1j*sigma_y)/2
I = Matrix([[1,0],
            [0,1]])
H = (omega*sigma_z + ux*sigma_x + uy*sigma_y)/2
rho = (tau*I + x*sigma_x + y*sigma_y + z*sigma_z)/2

rho_pt = (tau_pt*I + x_pt*sigma_x + y_pt*sigma_y + z_pt*sigma_z)/2
partie_droite = -1j*(H*rho - rho*H) + gamma*(sigma_plus*rho*sigma_moins - (sigma_moins*sigma_plus*rho + rho*sigma_moins*sigma_plus)/2)

equation = Eq(rho_pt,partie_droite)
solution = solve(equation,(tau_pt,x_pt,y_pt,z_pt))

print(solution)