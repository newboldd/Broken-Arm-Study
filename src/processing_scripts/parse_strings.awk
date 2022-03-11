BEGIN {
	go = 0;
	k = 0;
}

(NF==1 && go>1) {T[k] = $1; k++;}

/^MosaicRef/{go++;}
($1==0 && go==1) {T[k] = $1; k++; go++}
/AutoInline/{go=0;}

END {
	for (i = 0; i < k; i++) {
		if (T[i] == 0 && i > 0) exit;
		printf("%s\n",T[i]);
	}
}
