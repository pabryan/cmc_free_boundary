︠eb1e086b-a819-4768-9716-7566c2e00d96s︠
#### Some helper functions #####

# Get coefficients along the boundary
def boundary_coefficients(d):
    return {k: d[k].substitute({r: 0}) for k in d.keys()}

def boundary_vector(v):
    return v.substitute({r: 0})

# Get matrix version of 2-tensor
def matrix_form(T):
    return matrix([[T[(1,1)], T[(1,2)]], [T[(2,1)], T[(2,2)]]])

# Nice LaTeX form of connection coefficients
def cleanup_one_zero(text, num):
    ret = ""
    if num == 0:
        ret = ""
    elif num == 1:
        ret = text
    elif num == -1:
        ret = "-" + text
    else:
        ret = str(num) + text

    return ret

def connection_latex(k, gamma):
    lookup = {(1, 1) : ("T", "T"), (1, 2) : ("T", "N"), (2, 1) : ("N", "T"), (2, 2) : ("N", "N")}
    l = lookup[k]
    Ttext = cleanup_one_zero("T", gamma[k + (1,)])
    Ntext = cleanup_one_zero("N", gamma[k + (2,)])

    if Ttext == "" and Ntext == "":
        text = "0"
    elif Ttext == "":
        text = Ntext
    elif Ntext == "":
        text = Ttext
    else:
        text = Ttext + " + " + Ntext

    return r"$\nabla_{" + l[0] + r"}" + l[1] + r" = " + text + r"$"

# LaTeX display of variables
def latex_display(name, val):
    return r"$" + name + r" = " + latex(val) + r"$"
︡e4d49057-43ba-4b23-8968-c1be4b1eddb8︡{"done":true}︡
︠7b2d1551-6361-43ff-9791-3b94a5225304s︠
#### Base variables ###
s, r = var('s, r', domain='real')
︡3014d438-8aab-4815-81aa-4b98fbc6936a︡{"done":true}
︠d7eed68c-aea9-4ec8-bb81-a44bdfca5036︠
# Symbolic functions for general formulae
g1 = function('gamma1')(s)
g2 = function('gamma2')(s)
g3 = function('gamma3')(s)
gamma = vector([g1, g2, g3])
f = function('f')(s)

assume('g1(s)^2 + g2(s)^2 + g3(s)^3 = 1') # not sure if this does anything

# The boundary orthongality condition implies the inner normal is -gamma
# Our surface is constructed by following curves starting on gamma and over the inner normal
# where "over" means in the plane including the inner normal and orthogonal to gamma

# Computation is very slow using the general formulae and we seem to run out of memory!
gammap = gamma.diff(s)
V = gamma.cross_product(gammap)

# Computation is okay using
V = vector([0, 0, 1])

F = (1-r) * gamma + r^2 * f * V
coords = ((s, 0, 2*pi), (r, 0, 1/2))

︡8ffe2ee5-87b7-4e0e-851d-c1b2b8c068d3︡{"done":true}︡
︠a69a576b-3a0c-4b9a-bc03-dd2c23682c6ds︠
# Explicit functions for nice formulae and plots
R = 9/10
epsilon = 1/2

gamma = vector([R * cos(s/R), R * sin(s/R), sqrt(1-R^2)])
f = 2*cos(s/R)^3

gammap = gamma.diff(s)
V = gamma.cross_product(gammap)
F = (1-r) * gamma + r^2 * f * V
coords = ((s, 0, 2*pi*R), (r, 0, epsilon))
︡53355c93-16ad-4203-aa9f-066b84b8046b︡{"done":true}︡
︠00d7e8af-222f-40f3-a590-774bd16e4affs︠
#### Create our surface ####
M = ParametrizedSurface3D(F, coords, 'M')

︡452dc1c0-ef91-4c55-a13b-ec3db023d300︡{"done":true}︡
︠52c342cb-5236-426c-8f0f-abac5d103e9fs︠
#### Get fundamental quantities ####

# Coordinate frame
frame = M.natural_frame()
es = frame[1]
er = frame[2]

bdry_frame = boundary_coefficients(frame)
T = bdry_frame[1]
N = bdry_frame[2]

# Normal vector
n = M.normal_vector(normalized=True)
bdry_n = boundary_vector(n)

# Metric
g = M.first_fundamental_form_coefficients()
bdry_g = boundary_coefficients(g)
bdry_gmatrix = matrix_form(bdry_g)

# Connection
Gamma = M.connection_coefficients()
bdry_Gamma = boundary_coefficients(Gamma)

# Extrinsic curvature
h = M.second_fundamental_form_coefficients()
bdry_h = boundary_coefficients(h)
bdry_hmatrix = matrix_form(bdry_h)


︡37149d10-f821-4776-840b-6b8bfa911567︡{"done":true}︡
︠a56baa12-51d3-4208-a273-6490b2857db8s︠
#### Display fundamental quantities ####

print("Frame along boundary")
show(latex_display("T", T))
show(latex_display("N", N))

print("Normal along boundary")
show(latex_display("n", bdry_n))

print("Metric along boundary")
show(latex_display("g", bdry_gmatrix))

print("Connection along boundary")
show(connection_latex((1,1), bdry_Gamma))
show(connection_latex((1,2), bdry_Gamma))
show(connection_latex((2,1), bdry_Gamma))
show(connection_latex((2,2), bdry_Gamma))

print("Second fundamental form along boundary")
show(latex_display("h", bdry_hmatrix))
︡d72cf829-57d7-4446-8ea9-f260d47794ba︡{"stdout":"Frame along boundary\n"}︡{"html":"<div align='center'>$T = \\left(-\\sin\\left(\\frac{10}{9} \\, s\\right),\\,\\cos\\left(\\frac{10}{9} \\, s\\right),\\,0\\right) $</div>"}︡{"html":"<div align='center'>$N = \\left(-\\frac{9}{10} \\, \\cos\\left(\\frac{10}{9} \\, s\\right),\\,-\\frac{9}{10} \\, \\sin\\left(\\frac{10}{9} \\, s\\right),\\,-\\frac{1}{10} \\, \\sqrt{19}\\right) $</div>"}︡{"stdout":"Normal along boundary\n"}︡{"html":"<div align='center'>$n = \\left(-\\frac{1}{10} \\, \\sqrt{19} \\cos\\left(\\frac{10}{9} \\, s\\right),\\,-\\frac{1}{10} \\, \\sqrt{19} \\sin\\left(\\frac{10}{9} \\, s\\right),\\,\\frac{9}{10}\\right) $</div>"}︡{"stdout":"Metric along boundary\n"}︡{"html":"<div align='center'>$g = \\left(\\begin{array}{rr}\n1 &amp; 0 \\\\\n0 &amp; 1\n\\end{array}\\right) $</div>"}︡{"stdout":"Connection along boundary\n"}︡{"html":"<div align='center'>$\\nabla_{T}T = N$</div>"}︡{"html":"<div align='center'>$\\nabla_{T}N = -T$</div>"}︡{"html":"<div align='center'>$\\nabla_{N}T = -T$</div>"}︡{"html":"<div align='center'>$\\nabla_{N}N = 0$</div>"}︡{"stdout":"Second fundamental form along boundary\n"}︡{"html":"<div align='center'>$h = \\left(\\begin{array}{rr}\n\\frac{1}{9} \\, \\sqrt{19} &amp; 0 \\\\\n0 &amp; 4 \\, \\cos\\left(\\frac{10}{9} \\, s\\right)^{3}\n\\end{array}\\right) $</div>"}︡{"done":true}︡
︠36ec9645-a490-40ba-8542-7165fd895348s︠
#### Plot surface ####

p = M.plot(color='red') + sphere(opacity=0.3)
p.show()
︡4f171df9-03d8-4d63-8801-d3c74ad6989d︡{"file":{"filename":"cf54b9de-0088-4119-8ac0-4b16d37b478c.sage3d","uuid":"cf54b9de-0088-4119-8ac0-4b16d37b478c"}}︡{"done":true}︡
︠d94ad243-10e3-42d1-b28e-8974a69c3047︠

#### Check computations for variation of h(N, N) ####

## Lemma 5.4 ##
print("Lemma 5.4")

# kappa = -1
bool(bdry_Gamma[(1, 1, 1)] == 0)
bool(bdry_Gamma[(1, 1, 2)] == 1)

# T, N is an o/n basis
bool(bdry_gmatrix == matrix.identity(2))

# h is diagonalised by T, N
bool(bdry_h[(1, 2)] == 0)
bool(bdry_h[(2, 1)] == 0)

## Lemma 5.5 ##
print("Lemma 5.5")

Wr = -n.diff(r) # $Wr = -D_{\bar{N}} \nu^M$
WN = boundary_vector(Wr) # $WN = -D_N \nu^M$ along $\partial M$
Ws = -n.diff(s) # $Ws = -D_{\bar{T}} \nu^M$
WT = boundary_vector(Ws) # $WT = -D_T \nu^M$ along $\partial M$

# Equation (14)
bool(bdry_h[(2,2)].diff(s) == N * WN.diff(s))

# Commuting \bar{T}, \bar{N} in (14)
bool(er * Wr.diff(s) == er * Ws.diff(r))
bool(bdry_h[(2,2)].diff(s) == N * boundary_vector(Ws.diff(r))) # at r = 0

# Equation (18)
bool(boundary_vector(bdry_h[(2,2)].diff(s)) == boundary_vector(h[(1,2)].diff(r)))

︡482b5843-e345-4910-bff9-3bcc29fce9d2︡{"stdout":"Lemma 5.4\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"Lemma 5.5\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"stdout":"True\n"}︡{"done":true}︡
︠cbf3b365-b6a0-4470-afbf-af1a61b80948︠









