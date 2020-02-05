IFORT = ifort
MODS = -I/var/jobs/ppg06/LIB/modules
OPCOES = -warn -r8 -i8 -ip -O3 -axsse4.2 -u

MAIN = calc.lyap.mod.f90

all: invt ar

ar: invt
	ar vru calc.lyap.mod.a calc.lyap.mod.o

invt:
	$(IFORT) -c $(OPCOES) $(MAIN) $(lqr) $(MODS)