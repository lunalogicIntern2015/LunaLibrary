Bool decchk(string str, char &ch) {
	char c;
	Int j,k=0,m=0,n=str.length();
	static Int ip[10][8]={{0,1,5,8,9,4,2,7},{1,5,8,9,4,2,7,0},
		{2,7,0,1,5,8,9,4},{3,6,3,6,3,6,3,6},{4,2,7,0,1,5,8,9},
		{5,8,9,4,2,7,0,1},{6,3,6,3,6,3,6,3},{7,0,1,5,8,9,4,2},
		{8,9,4,2,7,0,1,5},{9,4,2,7,0,1,5,8}};
	static Int ij[10][10]={{0,1,2,3,4,5,6,7,8,9},{1,2,3,4,0,6,7,8,9,5},
		{2,3,4,0,1,7,8,9,5,6},{3,4,0,1,2,8,9,5,6,7},{4,0,1,2,3,9,5,6,7,8},
		{5,9,8,7,6,0,4,3,2,1},{6,5,9,8,7,1,0,4,3,2},{7,6,5,9,8,2,1,0,4,3},
		{8,7,6,5,9,3,2,1,0,4},{9,8,7,6,5,4,3,2,1,0}};
	for (j=0;j<n;j++) {
		c=str[j];
		if (c >= 48 && c <= 57)
			k=ij[k][ip[(c+2) % 10][7 & m++]];
	}
	for (j=0;j<10;j++)
		if (ij[k][ip[j][m & 7]] == 0) break;
	ch=char(j+48);
	return k==0;
}
