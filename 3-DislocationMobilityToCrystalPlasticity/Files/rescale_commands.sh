#!/bin/bash -x



#######################################################################
### define file variables and command line arguments
REFERENCE_STRUCTURE="bcc" # crystal structure
ELEMENT_NAME="Fe" # Periodic Table identifier of element
declare -a C=(294.55e3  144.472e3  195.343e3) # (C11, C12, C44), MPa
DRAG_COEFFICIENT=2.1576805328804272e-05 # Pa-s
H0=0.41169999999999873 # MPa
KAPPA_0=0.3101737291184037 # MPa
KAPPA_S=0.3101737291184037 # MPa
declare -a STRAIN_INCREMENT=(1000 1000 1000) # number of substeps
declare -a STRAIN_RATE=(0.0001 -0.0001 0.0001) # s
declare -a STRAIN=(1. 1. 1.) # number to repeat substeps




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




# reference structure selection
case $(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]') in

    "fcc")
        # elastic constants
        # $(printf '%g %g %g' 
        sed -i "2s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*%${C[0]} ${C[1]} ${C[2]}     / c11(c1), c12(c2), c44(c3) /%" "./fcc.sx"

        # drag coefficient
        sed -i "s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ bdrag /%$DRAG_COEFFICIENT     / bdrag /%" "./fcc.sx"

        # initial hardening rate
        # sed -i "s%[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ h0, tausi, taus0, xms, gamss0 /%$H0 3.7  30.80  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./fcc.sx"
        sed -i "14s%[[:print:]]*%$H0 $KAPPA_0  $KAPPA_S  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./fcc.sx"
        sed -i "34s%[[:print:]]*%$H0 $KAPPA_0  $KAPPA_S  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./fcc.sx"

        # initial strength
        sed -i "s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ initial slip system hardness (kappa0) /%$KAPPA_0     / initial slip system hardness (kappa0) /%" "./fcc.sx"
        ;;

    "bcc")
        # elastic constants
        # $(printf '%g %g %g' 
        sed -i "2s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*%${C[0]} ${C[1]} ${C[2]}     / c11(c1), c12(c2), c44(c3) /%" "./bcc.sx"

        # drag coefficient
        sed -i "s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ bdrag /%$DRAG_COEFFICIENT     / bdrag /%" "./bcc.sx"

        # initial hardening rate
        # sed -i "s%[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ h0, tausi, taus0, xms, gamss0 /%$H0 3.7  30.80  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./bcc.sx"
        sed -i "14s%[[:print:]]*%$H0 $KAPPA_0  $KAPPA_S  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./bcc.sx"
        sed -i "34s%[[:print:]]*%$H0 $KAPPA_0  $KAPPA_S  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./bcc.sx"
        sed -i "54s%[[:print:]]*%$H0 $KAPPA_0  $KAPPA_S  0.0-4  5.0e10     / h0, tausi, taus0, xms, gamss0 /%" "./bcc.sx"

        # initial strength
        sed -i "s%^[[:alpha:]]*.*[[:alpha:]]*e*-*[[:alpha:]]*[[:blank:]]*/ initial slip system hardness (kappa0) /%$KAPPA_0     / initial slip system hardness (kappa0) /%" "./bcc.sx"
        ;;
esac

# tension
mkdir "./tension"
cp *.* "./tension/"
cd "./tension"
sed -i "s%^incrmt = [[:alpha:]]*.*[[:alpha:]]*%incrmt = ${STRAIN_INCREMENT[0]}%" "./rve.single.inp"
sed -i "s%^rate   = [[:alpha:]]*.*[[:alpha:]]*%rate   = ${STRAIN_RATE[0]}%" "./rve.single.inp"
sed -i "s%^strain = [[:alpha:]]*.*[[:alpha:]]*%strain = ${STRAIN[0]}%" "./rve.single.inp"

sed -i "93s%[[:print:]]*%\*Boundary,user%" "./rve.single.inp"
sed -i "94s%\**[^[[:print:]]]*%%" "./rve.single.inp"
sed -i "96s%[[:print:]]*%\*Boundary%" "./rve.single.inp"
sed -i "97s%\**[^[[:print:]]]*%%" "./rve.single.inp"
sed -i "99s%[[:print:]]*%\*Boundary%" "./rve.single.inp"
sed -i "100s%\**[^[[:print:]]]*%%" "./rve.single.inp"
sed -i "102s%[[:print:]]*%\*Boundary%" "./rve.single.inp"
sed -i "103s%\**[^[[:print:]]]*%%" "./rve.single.inp"
sed -i "104s%\**[^[[:print:]]]*%%" "./rve.single.inp"
sed -i "105s%\**[^[[:print:]]]*%%" "./rve.single.inp"

PWD <<< pwd
sed -i "95s%^[[:print:]]*%     \&    /'$PWD'/%" "./umat_xtal.f"
sed -i "7828s%^      RATE = 0.0001D-00%      RATE = ${STRAIN_RATE[0]}D-00%" "./umat_xtal.f"

abaqus job=icme.cpfem.single.tension input=rve.single.inp user=umat_xtal.f cpus=2 double

cd "../"
mkdir "./compression"
cp *.* "./compression/"
cd "./compression"

# compression
sed -i "s%^incrmt = [[:alpha:]]*.*[[:alpha:]]*%incrmt = ${STRAIN_INCREMENT[1]}%" "./rve.single.inp"
sed -i "s%^rate   = [[:alpha:]]*.*[[:alpha:]]*%rate   = ${STRAIN_RATE[1]}%" "./rve.single.inp"
sed -i "s%^strain = [[:alpha:]]*.*[[:alpha:]]*%strain = ${STRAIN[1]}%" "./rve.single.inp"

PWD <<< pwd
sed -i "95s%^[[:print:]]*%     \&    /'$PWD'/%" "./umat_xtal.f"
sed -i "7828s%^      RATE = 0.0001D-00%      RATE = ${STRAIN_RATE[0]}D-00%" "./umat_xtal.f"

abaqus job=icme.cpfem.single.compression input=rve.single.inp user=umat_xtal.f cpus=2 double


cd "../"
mkdir "./shear"
cp *.* "./shear/"
cd "./shear"

# shear
sed -i "s%^incrmt = [[:alpha:]]*.*[[:alpha:]]*%incrmt = ${STRAIN_INCREMENT[2]}%" "./rve.single.inp"
sed -i "s%^rate   = [[:alpha:]]*.*[[:alpha:]]*%rate   = ${STRAIN_RATE[2]}%" "./rve.single.inp"
sed -i "s%^strain = [[:alpha:]]*.*[[:alpha:]]*%strain = ${STRAIN[2]}%" "./rve.single.inp"

sed -i "93s%[[:print:]]*%\*\*Boundary,user%" "./rve.single.inp"
sed -i "94s%\**grainone%\*\*grainone%" "./rve.single.inp"
sed -i "96s%[[:print:]]*%\*\*Boundary%" "./rve.single.inp"
sed -i "97s%\**grainone%\*\*grainone%" "./rve.single.inp"
sed -i "99s%[[:print:]]*%\*\*Boundary%" "./rve.single.inp"
sed -i "100s%\**grainone%\*\*grainone%" "./rve.single.inp"
sed -i "102s%[[:print:]]*%\*\*Boundary%" "./rve.single.inp"
sed -i "103s%\**grainone%\*\*grainone%" "./rve.single.inp"
sed -i "104s%\**grainone%\*\*grainone%" "./rve.single.inp"
sed -i "105s%\**grainone%\*\*grainone%" "./rve.single.inp"

sed -i "110s%[[:print:]]*%\*Boundary,user%" "./rve.single.inp"
sed -i "111s%\**grainone%grainone%" "./rve.single.inp"
sed -i "113s%[[:print:]]*%\*Boundary,user%" "./rve.single.inp"
sed -i "114s%\**grainone%grainone%" "./rve.single.inp"
sed -i "116s%[[:print:]]*%\*Boundary,user%" "./rve.single.inp"
sed -i "117s%\**grainone%grainone%" "./rve.single.inp"

PWD <<< pwd
sed -i "95s%^[[:print:]]*%     \&    /'$PWD'/%" "./umat_xtal.f"
sed -i "7828s%^      RATE = 0.0001D-00%      RATE = ${STRAIN_RATE[0]}D-00%" "./umat_xtal.f"

abaqus job=icme.cpfem.single.shear input=rve.single.inp user=umat_xtal.f cpus=2 double




# end of file