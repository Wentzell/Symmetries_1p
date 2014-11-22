#!/bin/bash
Pi=3.14159265
PAR[1]="0.5" #G_L	PAR[1]="0.5" #G_L			
PAR[2]="10.0" #D      PAR[2]="1e+3" #D		    
PAR[3]="3.0" #U		PAR[3]="3.14159" #U  
PAR[4]="0.0" #E         PAR[4]="0.1" #E			    
PAR[5]="0.1" #PHI/PI    PAR[5]="1.0" #PHI/PI	    
PAR[6]="0.0" #B         PAR[6]="0.0" #B			    
PAR[7]="10.0" #beta

PAR[8]="10:00" #Required Time
PAR[9]="1e-10" #Error
ans="7"

fname=GL"${PAR[1]}"_D"${PAR[2]}"_U"${PAR[3]}"_E"${PAR[4]}"_Phi"${PAR[5]}"Pi_B"${PAR[6]}"_Beta"${PAR[7]}" # CREATE PART OF FILENAME

rm dat*.sh #dat${name[$i]}$fname.sh
echo "rm callgrind.out.*" >> dat_$fname.sh
echo "time nice -n 19 valgrind --tool=callgrind ./bin/run " >> dat_$fname.sh
#\"${PAR[1]}\" \"${PAR[2]}\" \"${PAR[3]}\" \"${PAR[4]}\" \"${PAR[5]}\" \"${PAR[6]}\" \"${PAR[7]}\" \"$fname.dat\" \"${PAR[9]}\" 
echo "callgrind_annotate --auto=yes callgrind.out.* >> annot_src.cpp" >> dat_$fname.sh
echo "kcachegrind callgrind.out.* &" >> dat_$fname.sh
chmod a+x dat_$fname.sh
./dat_$fname.sh
