awk '{OFS="\t"}{ 
  if(index($1,"@")==1) { 
    print;
  } else {
    match($6,"^([0-9]+)S",fiveclip); 
    if(length(fiveclip)==0){
       fiveclip[1]=0;
    }
    match($6,"([0-9]+)S$",threeclip);
    if(length(threeclip)==0){
       threeclip[1]=0;
    }
    # parse intron junction motif tag (jM) from STAR SAM file format output
    # score >= 20 for those reads that span annotated junction
    split($18,a,","); 
    isAnnoSplice = 0; 
    for(i=2;i<=length(a);i++){ 
       if(a[i]>=20) isAnnoSplice=1;
    } 
    if(isAnnoSplice>0 &&  fiveclip[1]<5 && threeclip[1]<5) print ;
  }
}'
