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
︡d76043ec-01bb-4bd2-bba0-904faafda4dc︡{"done":true}︡
︠7b2d1551-6361-43ff-9791-3b94a5225304s︠
#### Base variables ###
s, r = var('s, r', domain='real')
︡81adbcd5-b5a1-469a-83cd-986d7d12345a︡{"done":true}︡
︠d7eed68c-aea9-4ec8-bb81-a44bdfca5036s︠
# Symbolic functions for general formulae
g1 = function('gamma1')(s)
g2 = function('gamma2')(s)
g3 = function('gamma3')(s)
gamma = vector([g1, g2, g3])
f = function('f')(s)

assume('g1(s)^2 + g2(s)^2 + g3(s)^3 = 1')
︡eaa79a86-bd9c-4494-a918-51329f06ba90︡{"done":true}︡
︠a69a576b-3a0c-4b9a-bc03-dd2c23682c6ds︠
# Explicit functions for nice formulae and plots
gamma = vector([cos(s), sin(s), 0])
f = 2*cos(s)^3
︡92357b37-6ee6-47b1-a6ae-71c879dcbaec︡{"done":true}︡
︠00d7e8af-222f-40f3-a590-774bd16e4affs︠
#### Create our surface ####

F = (1-r) * gamma + vector([0, 0, f * r^2])
coords = ((s, 0, 2*pi), (r, 0, 1/2))
M = ParametrizedSurface3D(F, coords, 'M')

︡9fe5061b-a7bb-46c1-a484-fd0a26f754a7︡{"done":true}︡
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


︡da7d721f-1fb0-4b64-a96b-3142d4b8c514︡{"done":true}︡
︠a56baa12-51d3-4208-a273-6490b2857db8s︠
#### Plot surface and display fundamental quantities ####

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

# Plot
#p = M.plot(color='red') + sphere(opacity=0.3)
#p.show()
︡6b402bf3-f039-4a6d-a786-d45807888769︡{"stdout":"Frame along boundary\n"}︡{"html":"<div align='center'>$T = \\left(-\\sin\\left(s\\right),\\,\\cos\\left(s\\right),\\,0\\right) $</div>"}︡{"html":"<div align='center'>$N = \\left(-\\cos\\left(s\\right),\\,-\\sin\\left(s\\right),\\,0\\right) $</div>"}︡{"stdout":"Normal along boundary\n"}︡{"html":"<div align='center'>$n = \\left(0,\\,0,\\,1\\right) $</div>"}︡{"stdout":"Metric along boundary\n"}︡{"html":"<div align='center'>$g = \\left(\\begin{array}{rr}\n1 &amp; 0 \\\\\n0 &amp; 1\n\\end{array}\\right) $</div>"}︡{"stdout":"Connection along boundary\n"}︡{"html":"<div align='center'>$\\nabla_{T}T = N$</div>"}︡{"html":"<div align='center'>$\\nabla_{T}N = -T$</div>"}︡{"html":"<div align='center'>$\\nabla_{N}T = -T$</div>"}︡{"html":"<div align='center'>$\\nabla_{N}N = 0$</div>"}︡{"stdout":"Second fundamental form along boundary\n"}︡{"html":"<div align='center'>$h = \\left(\\begin{array}{rr}\n0 &amp; 0 \\\\\n0 &amp; 4 \\, \\cos\\left(s\\right)^{3}\n\\end{array}\\right) $</div>"}︡{"done":true}︡
︠d94ad243-10e3-42d1-b28e-8974a69c3047s︠

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









