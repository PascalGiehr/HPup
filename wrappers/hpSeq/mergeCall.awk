BEGIN{
 FS="\t";
 OFS="\t";
}

{
 key=$1"@"$2;
}

key!=curKey{
 if(NR>1) {
  print chr,pos,r1,r2;
 } 
 curKey=key;
 chr=$1;
 pos=$2;
 r1="";
 r2="";
}

{ 
r1=r1$3;
r2=r2$4;
} 

END{
 print chr,pos,r1,r2;
}
