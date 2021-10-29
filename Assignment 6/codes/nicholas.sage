
N = 84364443735725034864402554533826279174703893439763343343863260342756678609216895093779263028809246505955647572176682669445270008816481771701417554768871285020442403001649254405058303439906229201909599348669565697534331652019516409514800265887388539283381053937433496994442146419682027649079704982600857517093
C = 23701787746829110396789094907319830305538180376427283226295906585301889543996533410539381779684366880970896279018807100530176651625086988655210858554133345906272561027798171440923147960165094891980452757852685707020289384698322665347609905744582248157246932007978339129630067022987966706955482598869800151693
e = 5

#set of hexadecimal values obtained while trying different exit commands
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

print(f"Message : '{message}'")


ZmodN = IntegerModRing(N)
P.<R> = PolynomialRing(ZmodN)
h = 3
X = int((N**(1/9).n()))
k = 5
mat_dim = h*k


for l in range(1,200):
    pol = ((int_message * (2**l)) + R)**e - C
    pol = pol.change_ring(ZZ)
    x = pol.parent().gen()

    m = matrix(ZZ, mat_dim, mat_dim)

    for i in range(0,mat_dim):
        v = (i//k)
        u = i - k*v
        for j in range(0,mat_dim):
            coeff = (x**u) * (pol**v)
            e_ij = coeff[j] * (N**(h-1-v))
            m[i,j] =  (e_ij) * (X**j)

    m = m.LLL()

    m = m.transpose()
    for i in range(m.nrows()):
        m[i] = m[i]/(X**i)
    m = m.transpose()

    poly = 0

    P.<pol1> = PolynomialRing(ZZ)
    for j in range(m.nrows()):
        poly += m[0,j]*(pol1**j)


    #find the roots
    root = poly.roots()

    #if root exists
    if root:
        root = root[0][0]

        bin_root = ""
        temp = root

        #convert it into binary
        while(temp!=1):
            rem = temp % 2
            bin_root = str(rem) + bin_root
            temp = temp//2
        bin_root = "01"+bin_root

        #decode password by converting binary root using ASCII char code taking 8 bits at a time
        password = ""
        for i in range(0,len(bin_root),8):
            password += chr(int(bin_root[i:i+8],2))

        print(f"l = {l} \nroot : {root}  \nbinary root : {bin_root}")
        print(f"PASSWORD ::: {password}")
        break

print("DONE")











