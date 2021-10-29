N = 84364443735725034864402554533826279174703893439763343343863260342756678609216895093779263028809246505955647572176682669445270008816481771701417554768871285020442403001649254405058303439906229201909599348669565697534331652019516409514800265887388539283381053937433496994442146419682027649079704982600857517093
C = 23701787746829110396789094907319830305538180376427283226295906585301889543996533410539381779684366880970896279018807100530176651625086988655210858554133345906272561027798171440923147960165094891980452757852685707020289384698322665347609905744582248157246932007978339129630067022987966706955482598869800151693
e = 5

t=["59 6f 75 20 73 65 65 20","61 20 47 6f 6c 64 2d 42","75 67 20 69 6e 20 6f 6e","65 20 63 6f 72 6e 65 72","2e 20 49 74 20 69 73 20",
   "74 68 65 20 6b 65 79 20","74 6f 20 61 20 74 72 65","61 73 75 72 65 20 66 6f","75 6e 64 20 62 79"]

message=""
bin_message=""

for y in t:
    for x in y.split(" "):
        binary_byte = "{:0>8b}".format(int(x,16))
        bin_message +=binary_byte
        message+=chr(int(binary_byte,2))

#add one space at the end
message +=" "
bin_message += "{:0>8b}".format(ord(" "))
int_message = int(bin_message,2)

print(f"Message : {message}")

############### RSA ################

def coppersmith_howgrave_univariate(pol,dd, modulus, beta, mm, tt, XX):

    nn = dd * mm + tt
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)

    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()

    # test roots
    roots = []
    for root in potential_roots:
        if root[0].is_integer():
            result = polZ(ZZ(root[0]))
            if gcd(modulus, result) >= modulus^beta:
                roots.append(ZZ(root[0]))

    return roots


ZmodN = Zmod(N)
P.<x> = PolynomialRing(ZmodN)

for l in range(1,200):
    pol = ((int_message<<l) + x)^e - C
    dd = pol.degree()

    beta = 1
    epsilon = beta / 7
    mm = ceil(beta**2 / (dd * epsilon))
    tt = floor(dd * mm * ((1/beta) - 1))
    XX = ceil(N**((beta**2/dd) - epsilon))

    root = coppersmith_howgrave_univariate(pol,dd, N, beta, mm, tt, XX)
    if root:
        bin_root = "0" + root[0].binary()
        password = ""
        for i in range(0,len(bin_root),8):
            password += chr(int(bin_root[i:i+8],2))

        print(f"l = {l} \nroot : {root[0]}  \nbinary root : {bin_root}")
        print(f"PASSWORD ::: {password}")
        break

print("DONE")


