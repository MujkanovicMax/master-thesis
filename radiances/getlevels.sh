ATM="/home/m/Mujkanovic.Max/ma/radiances/stdatm/afglus_mgl_full.dat"
ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
NLAY=$( cat $ATM | tail -n +4  | wc -l )
ZLEVno0=$( echo $ZLEV | sed "s/^[^ ]* /-999 /" )

echo $ZLEVno0
