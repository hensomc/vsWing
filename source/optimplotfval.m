%!PS-Adobe-2.0 
%%Creator: MATLAB, The MathWorks, Inc. Version 8.3.0.532 (R2014a). Operating System: Microsoft Windows 7 Home Premium .
%%Title: .\optimplotfval.m
%%CreationDate: 03/10/2017  16:39:14
%%DocumentNeededFonts: Helvetica
%%DocumentProcessColors: Cyan Magenta Yellow Black
%%Pages: (atend)
%%BoundingBox: (atend)
%%EndComments

%%BeginProlog
% MathWorks dictionary
/MathWorks 160 dict begin
% definition operators
/bdef {bind def} bind def
/ldef {load def} bind def
/xdef {exch def} bdef
/xstore {exch store} bdef
% operator abbreviations
/c  /clip ldef
/cc /concat ldef
/cp /closepath ldef
/gr /grestore ldef
/gs /gsave ldef
/mt /moveto ldef
/np /newpath ldef
/cm /currentmatrix ldef
/sm /setmatrix ldef
/rm /rmoveto ldef
/rl /rlineto ldef
/s {show newpath} bdef
/sc {setcmykcolor} bdef
/sr /setrgbcolor ldef
/sg /setgray ldef
/w /setlinewidth ldef
/j /setlinejoin ldef
/cap /setlinecap ldef
/rc {rectclip} bdef
/rf {rectfill} bdef
% page state control
/pgsv () def
/bpage {/pgsv save def} bdef
/epage {pgsv restore} bdef
/bplot /gsave ldef
/eplot {stroke grestore} bdef
% orientation switch
/portraitMode 0 def /landscapeMode 1 def /rotateMode 2 def
% coordinate system mappings
/dpi2point 0 def
% font control
/FontSize 0 def
/FMS {/FontSize xstore findfont [FontSize 0 0 FontSize neg 0 0]
  makefont setfont} bdef
/ISOLatin1Encoding where {pop /WindowsLatin1Encoding 256 array bdef
ISOLatin1Encoding WindowsLatin1Encoding copy pop
/.notdef/.notdef/quotesinglbase/florin/quotedblbase/ellipsis/dagger
/daggerdbl/circumflex/perthousand/Scaron/guilsinglleft/OE/.notdef/.notdef
/.notdef/.notdef/quoteleft/quoteright/quotedblleft/quotedblright/bullet
/endash/emdash/tilde/trademark/scaron/guilsinglright/oe/.notdef/.notdef
/Ydieresis WindowsLatin1Encoding 128 32 getinterval astore pop}
{/WindowsLatin1Encoding StandardEncoding bdef} ifelse
/reencode {exch dup where {pop load} {pop StandardEncoding} ifelse
  exch dup 3 1 roll findfont dup length dict begin
  { 1 index /FID ne {def}{pop pop} ifelse } forall
  /Encoding exch def currentdict end definefont pop} bdef
/isroman {findfont /CharStrings get /Agrave known} bdef
/FMSR {3 1 roll 1 index dup isroman {reencode} {pop pop} ifelse
  exch FMS} bdef
/csm {1 dpi2point div -1 dpi2point div scale neg translate
 dup landscapeMode eq {pop -90 rotate}
  {rotateMode eq {90 rotate} if} ifelse} bdef
% line types: solid, dotted, dashed, dotdash
/SO { [] 0 setdash } bdef
/DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef
/DA { [6 dpi2point mul] 0 setdash } bdef
/DD { [.5 dpi2point mul 4 dpi2point mul 6 dpi2point mul 4
  dpi2point mul] 0 setdash } bdef
% macros for lines and objects
/L {lineto stroke} bdef
/MP {3 1 roll moveto 1 sub {rlineto} repeat} bdef
/AP {{rlineto} repeat} bdef
/PDlw -1 def
/W {/PDlw currentlinewidth def setlinewidth} def
/PP {closepath eofill} bdef
/DP {closepath stroke} bdef
/MR {4 -2 roll moveto dup  0 exch rlineto exch 0 rlineto
  neg 0 exch rlineto closepath} bdef
/FR {MR stroke} bdef
/PR {MR fill} bdef
/L1i {{currentfile picstr readhexstring pop} image} bdef
/tMatrix matrix def
/MakeOval {newpath tMatrix currentmatrix pop translate scale
0 0 1 0 360 arc tMatrix setmatrix} bdef
/FO {MakeOval stroke} bdef
/PO {MakeOval fill} bdef
/PD {currentlinewidth 2 div 0 360 arc fill
   PDlw -1 eq not {PDlw w /PDlw -1 def} if} def
/FA {newpath tMatrix currentmatrix pop translate scale
  0 0 1 5 -2 roll arc tMatrix setmatrix stroke} bdef
/PA {newpath tMatrix currentmatrix pop	translate 0 0 moveto scale
  0 0 1 5 -2 roll arc closepath tMatrix setmatrix fill} bdef
/FAn {newpath tMatrix currentmatrix pop translate scale
  0 0 1 5 -2 roll arcn tMatrix setmatrix stroke} bdef
/PAn {newpath tMatrix currentmatrix pop translate 0 0 moveto scale
  0 0 1 5 -2 roll arcn closepath tMatrix setmatrix fill} bdef
/vradius 0 def /hradius 0 def /lry 0 def
/lrx 0 def /uly 0 def /ulx 0 def /rad 0 def
/MRR {/vradius xdef /hradius xdef /lry xdef /lrx xdef /uly xdef
  /ulx xdef newpath tMatrix currentmatrix pop ulx hradius add uly
  vradius add translate hradius vradius scale 0 0 1 180 270 arc 
  tMatrix setmatrix lrx hradius sub uly vradius add translate
  hradius vradius scale 0 0 1 270 360 arc tMatrix setmatrix
  lrx hradius sub lry vradius sub translate hradius vradius scale
  0 0 1 0 90 arc tMatrix setmatrix ulx hradius add lry vradius sub
  translate hradius vradius scale 0 0 1 90 180 arc tMatrix setmatrix
  closepath} bdef
/FRR {MRR stroke } bdef
/PRR {MRR fill } bdef
/MlrRR {/lry xdef /lrx xdef /uly xdef /ulx xdef /rad lry uly sub 2 div def
  newpath tMatrix currentmatrix pop ulx rad add uly rad add translate
  rad rad scale 0 0 1 90 270 arc tMatrix setmatrix lrx rad sub lry rad
  sub translate rad rad scale 0 0 1 270 90 arc tMatrix setmatrix
  closepath} bdef
/FlrRR {MlrRR stroke } bdef
/PlrRR {MlrRR fill } bdef
/MtbRR {/lry xdef /lrx xdef /uly xdef /ulx xdef /rad lrx ulx sub 2 div def
  newpath tMatrix currentmatrix pop ulx rad add uly rad add translate
  rad rad scale 0 0 1 180 360 arc tMatrix setmatrix lrx rad sub lry rad
  sub translate rad rad scale 0 0 1 0 180 arc tMatrix setmatrix
  closepath} bdef
/FtbRR {MtbRR stroke } bdef
/PtbRR {MtbRR fill } bdef
/stri 6 array def /dtri 6 array def
/smat 6 array def /dmat 6 array def
/tmat1 6 array def /tmat2 6 array def /dif 3 array def
/asub {/ind2 exch def /ind1 exch def dup dup
  ind1 get exch ind2 get sub exch } bdef
/tri_to_matrix {
  2 0 asub 3 1 asub 4 0 asub 5 1 asub
  dup 0 get exch 1 get 7 -1 roll astore } bdef
/compute_transform {
  dmat dtri tri_to_matrix tmat1 invertmatrix 
  smat stri tri_to_matrix tmat2 concatmatrix } bdef
/ds {stri astore pop} bdef
/dt {dtri astore pop} bdef
/db {2 copy /cols xdef /rows xdef mul dup string
  currentfile exch readhexstring pop
  /bmap xdef pop pop} bdef
/it {gs np dtri aload pop moveto lineto lineto cp c
  cols rows 8 compute_transform 
  {bmap} image gr}bdef
/il {newpath moveto lineto stroke}bdef
currentdict end def
%%EndProlog

%%BeginSetup
MathWorks begin

0 cap

end
%%EndSetup

%%Page: 1 1
%%BeginPageSetup
%%PageBoundingBox:    66   212   567   583
MathWorks begin
bpage
%%EndPageSetup

%%BeginObject: obj1
bplot

/dpi2point 12 def
portraitMode 0204 7344 csm

  589   340  6017  4460 MR c np
85 dict begin %Colortable dictionary
/c0 { 0.000000 0.000000 0.000000 sr} bdef
/c1 { 1.000000 1.000000 1.000000 sr} bdef
/c2 { 0.900000 0.000000 0.000000 sr} bdef
/c3 { 0.000000 0.820000 0.000000 sr} bdef
/c4 { 0.000000 0.000000 0.800000 sr} bdef
/c5 { 0.910000 0.820000 0.320000 sr} bdef
/c6 { 1.000000 0.260000 0.820000 sr} bdef
/c7 { 0.000000 0.820000 0.820000 sr} bdef
c0
1 j
1 sg
   0    0 6913 5185 PR
6 w
0 2335 -2325 -1070 0 -2335 6254 3793 4 MP
PP
2325 1070 0 2335 -2325 -1070 0 -2335 6254 3793 5 MP stroke
2325 1070 3031 -821 -2326 -1069 3224 4613 4 MP
PP
-3030 820 2325 1070 3031 -821 -2326 -1069 3224 4613 5 MP stroke
0 2335 3031 -821 0 -2335 898 3544 4 MP
PP
-3031 821 0 2335 3031 -821 0 -2335 898 3544 5 MP stroke
4 w
DO
0 sg
3224 4613 mt  898 3544 L
 898 3544 mt  898 1209 L
4234 4340 mt 1909 3270 L
1909 3270 mt 1909  936 L
5244 4066 mt 2919 2997 L
2919 2997 mt 2919  662 L
6254 3793 mt 3929 2723 L
3929 2723 mt 3929  388 L
3224 4613 mt 6254 3793 L
6254 3793 mt 6254 1458 L
2643 4346 mt 5673 3525 L
5673 3525 mt 5673 1191 L
2061 4079 mt 5092 3258 L
5092 3258 mt 5092  923 L
1480 3811 mt 4510 2991 L
4510 2991 mt 4510  656 L
 898 3544 mt 3929 2723 L
3929 2723 mt 3929  388 L
 898 3544 mt 3929 2723 L
3929 2723 mt 6254 3793 L
 898 3155 mt 3929 2334 L
3929 2334 mt 6254 3404 L
 898 2766 mt 3929 1945 L
3929 1945 mt 6254 3015 L
 898 2377 mt 3929 1556 L
3929 1556 mt 6254 2625 L
 898 1987 mt 3929 1167 L
3929 1167 mt 6254 2236 L
 898 1598 mt 3929  778 L
3929  778 mt 6254 1847 L
 898 1209 mt 3929  388 L
3929  388 mt 6254 1458 L
SO
6 w
3224 4613 mt 6254 3793 L
3224 4613 mt  898 3544 L
 898 3544 mt  898 1209 L
3224 4613 mt 3295 4646 L
%%IncludeResource: font Helvetica
/Helvetica /WindowsLatin1Encoding 120 FMSR

3327 4772 mt 
(0) s
4234 4340 mt 4305 4373 L
4337 4498 mt 
(50) s
5244 4066 mt 5316 4099 L
5347 4225 mt 
(100) s
6254 3793 mt 6326 3826 L
6357 3951 mt 
(150) s
3224 4613 mt 3148 4634 L
3049 4754 mt 
(0) s
2643 4346 mt 2567 4367 L
2400 4487 mt 
(50) s
2061 4079 mt 1985 4099 L
1752 4219 mt 
(100) s
1480 3811 mt 1404 3832 L
1171 3952 mt 
(150) s
 898 3544 mt  823 3565 L
 589 3685 mt 
(200) s
 898 3544 mt  827 3511 L
 660 3541 mt 
(-6) s
 898 3155 mt  827 3122 L
 660 3152 mt 
(-5) s
 898 2766 mt  827 2733 L
 660 2763 mt 
(-4) s
 898 2377 mt  827 2344 L
 660 2373 mt 
(-3) s
 898 1987 mt  827 1955 L
 660 1984 mt 
(-2) s
 898 1598 mt  827 1566 L
 660 1595 mt 
(-1) s
 898 1209 mt  827 1176 L
 730 1206 mt 
(0) s
gs 898 388 5357 4226 MR c np
/c8 { 0.500000 0.000000 0.000000 sr} bdef
c8
-83 276 -144 -288 4679 1885 3 MP
PP
0 sg
4679 1885 mt 4535 1597 L
4535 1597 mt 4452 1873 L
c8
-134 -338 -93 326 4679 1885 3 MP
PP
0 sg
4679 1885 mt 4586 2211 L
4586 2211 mt 4452 1873 L
/c9 { 1.000000 0.062500 0.000000 sr} bdef
c9
-83 260 -134 -338 4586 2211 3 MP
PP
0 sg
4586 2211 mt 4452 1873 L
4452 1873 mt 4369 2133 L
c9
-125 -386 -92 308 4586 2211 3 MP
PP
0 sg
4586 2211 mt 4494 2519 L
4494 2519 mt 4369 2133 L
c8
-74 225 -144 -287 4535 1597 3 MP
PP
0 sg
4535 1597 mt 4391 1310 L
4391 1310 mt 4317 1535 L
c8
-135 -338 -83 276 4535 1597 3 MP
PP
0 sg
4535 1597 mt 4452 1873 L
4452 1873 mt 4317 1535 L
/c10 { 1.000000 0.562500 0.000000 sr} bdef
c10
-116 -431 -92 290 4494 2519 3 MP
PP
0 sg
4494 2519 mt 4402 2809 L
4402 2809 mt 4286 2378 L
c10
-83 245 -125 -386 4494 2519 3 MP
PP
0 sg
4494 2519 mt 4369 2133 L
4369 2133 mt 4286 2378 L
/c11 { 0.937500 0.000000 0.000000 sr} bdef
c11
-126 -386 -83 260 4452 1873 3 MP
PP
0 sg
4452 1873 mt 4369 2133 L
4369 2133 mt 4243 1747 L
c11
-74 212 -135 -338 4452 1873 3 MP
PP
0 sg
4452 1873 mt 4317 1535 L
4317 1535 mt 4243 1747 L
/c12 { 0.937500 1.000000 0.062500 sr} bdef
c12
-107 -473 -93 272 4402 2809 3 MP
PP
0 sg
4402 2809 mt 4309 3081 L
4309 3081 mt 4202 2608 L
c12
-84 230 -116 -431 4402 2809 3 MP
PP
0 sg
4402 2809 mt 4286 2378 L
4286 2378 mt 4202 2608 L
c8
-65 174 -144 -287 4391 1310 3 MP
PP
0 sg
4391 1310 mt 4247 1023 L
4247 1023 mt 4182 1197 L
c8
-135 -338 -74 225 4391 1310 3 MP
PP
0 sg
4391 1310 mt 4317 1535 L
4317 1535 mt 4182 1197 L
/c13 { 1.000000 0.375000 0.000000 sr} bdef
c13
-116 -431 -83 245 4369 2133 3 MP
PP
0 sg
4369 2133 mt 4286 2378 L
4286 2378 mt 4170 1947 L
c13
-73 200 -126 -386 4369 2133 3 MP
PP
0 sg
4369 2133 mt 4243 1747 L
4243 1747 mt 4170 1947 L
/c14 { 0.875000 0.000000 0.000000 sr} bdef
c14
-64 164 -135 -338 4317 1535 3 MP
PP
0 sg
4317 1535 mt 4182 1197 L
4182 1197 mt 4118 1361 L
c14
-125 -386 -74 212 4317 1535 3 MP
PP
0 sg
4317 1535 mt 4243 1747 L
4243 1747 mt 4118 1361 L
/c15 { 0.437500 1.000000 0.562500 sr} bdef
c15
-97 -512 -93 253 4309 3081 3 MP
PP
0 sg
4309 3081 mt 4216 3334 L
4216 3334 mt 4119 2822 L
c15
-83 214 -107 -473 4309 3081 3 MP
PP
0 sg
4309 3081 mt 4202 2608 L
4202 2608 mt 4119 2822 L
/c16 { 1.000000 0.812500 0.000000 sr} bdef
c16
-106 -473 -84 230 4286 2378 3 MP
PP
0 sg
4286 2378 mt 4202 2608 L
4202 2608 mt 4096 2135 L
c16
-74 188 -116 -431 4286 2378 3 MP
PP
0 sg
4286 2378 mt 4170 1947 L
4170 1947 mt 4096 2135 L
c8
-54 123 -145 -287 4247 1023 3 MP
PP
0 sg
4247 1023 mt 4102  736 L
4102  736 mt 4048  859 L
c8
-134 -338 -65 174 4247 1023 3 MP
PP
0 sg
4247 1023 mt 4182 1197 L
4182 1197 mt 4048  859 L
/c17 { 1.000000 0.250000 0.000000 sr} bdef
c17
-64 155 -125 -386 4243 1747 3 MP
PP
0 sg
4243 1747 mt 4118 1361 L
4118 1361 mt 4054 1516 L
c17
-116 -431 -73 200 4243 1747 3 MP
PP
0 sg
4243 1747 mt 4170 1947 L
4170 1947 mt 4054 1516 L
/c18 { 0.000000 1.000000 1.000000 sr} bdef
c18
-83 197 -97 -512 4216 3334 3 MP
PP
0 sg
4216 3334 mt 4119 2822 L
4119 2822 mt 4036 3019 L
c18
-87 -548 -93 233 4216 3334 3 MP
PP
0 sg
4216 3334 mt 4123 3567 L
4123 3567 mt 4036 3019 L
/c19 { 0.750000 1.000000 0.250000 sr} bdef
c19
-97 -512 -83 214 4202 2608 3 MP
PP
0 sg
4202 2608 mt 4119 2822 L
4119 2822 mt 4022 2310 L
c19
-74 175 -106 -473 4202 2608 3 MP
PP
0 sg
4202 2608 mt 4096 2135 L
4096 2135 mt 4022 2310 L
/c20 { 0.750000 0.000000 0.000000 sr} bdef
c20
-125 -386 -64 164 4182 1197 3 MP
PP
0 sg
4182 1197 mt 4118 1361 L
4118 1361 mt 3993  975 L
c20
-55 116 -134 -338 4182 1197 3 MP
PP
0 sg
4182 1197 mt 4048  859 L
4048  859 mt 3993  975 L
c10
-65 145 -116 -431 4170 1947 3 MP
PP
0 sg
4170 1947 mt 4054 1516 L
4054 1516 mt 3989 1661 L
c10
-107 -474 -74 188 4170 1947 3 MP
PP
0 sg
4170 1947 mt 4096 2135 L
4096 2135 mt 3989 1661 L
/c21 { 0.000000 0.625000 1.000000 sr} bdef
c21
-78 -581 -93 209 4123 3567 3 MP
PP
0 sg
4123 3567 mt 4030 3776 L
4030 3776 mt 3952 3195 L
c21
-84 176 -87 -548 4123 3567 3 MP
PP
0 sg
4123 3567 mt 4036 3019 L
4036 3019 mt 3952 3195 L
c15
-74 160 -97 -512 4119 2822 3 MP
PP
0 sg
4119 2822 mt 4022 2310 L
4022 2310 mt 3948 2470 L
c15
-88 -549 -83 197 4119 2822 3 MP
PP
0 sg
4119 2822 mt 4036 3019 L
4036 3019 mt 3948 2470 L
c9
-55 110 -125 -386 4118 1361 3 MP
PP
0 sg
4118 1361 mt 3993  975 L
3993  975 mt 3938 1085 L
c9
-116 -431 -64 155 4118 1361 3 MP
PP
0 sg
4118 1361 mt 4054 1516 L
4054 1516 mt 3938 1085 L
c8
-45 72 -144 -287 4102 736 3 MP
PP
0 sg
4102  736 mt 3958  449 L
3958  449 mt 3913  521 L
c8
-135 -338 -54 123 4102 736 3 MP
PP
0 sg
4102  736 mt 4048  859 L
4048  859 mt 3913  521 L
/c22 { 1.000000 0.875000 0.000000 sr} bdef
c22
-97 -513 -74 175 4096 2135 3 MP
PP
0 sg
4096 2135 mt 4022 2310 L
4022 2310 mt 3925 1797 L
c22
-64 136 -107 -474 4096 2135 3 MP
PP
0 sg
4096 2135 mt 3989 1661 L
3989 1661 mt 3925 1797 L
/c23 { 1.000000 0.312500 0.000000 sr} bdef
c23
-107 -473 -65 145 4054 1516 3 MP
PP
0 sg
4054 1516 mt 3989 1661 L
3989 1661 mt 3882 1188 L
c23
-56 103 -116 -431 4054 1516 3 MP
PP
0 sg
4054 1516 mt 3938 1085 L
3938 1085 mt 3882 1188 L
/c24 { 0.687500 0.000000 0.000000 sr} bdef
c24
-126 -386 -55 116 4048 859 3 MP
PP
0 sg
4048  859 mt 3993  975 L
3993  975 mt 3867  589 L
c24
-46 68 -135 -338 4048 859 3 MP
PP
0 sg
4048  859 mt 3913  521 L
3913  521 mt 3867  589 L
/c25 { 0.062500 1.000000 0.937500 sr} bdef
c25
-74 145 -88 -549 4036 3019 3 MP
PP
0 sg
4036 3019 mt 3948 2470 L
3948 2470 mt 3874 2615 L
c25
-78 -580 -84 176 4036 3019 3 MP
PP
0 sg
4036 3019 mt 3952 3195 L
3952 3195 mt 3874 2615 L
/c26 { 0.000000 0.250000 1.000000 sr} bdef
c26
-84 154 -78 -581 4030 3776 3 MP
PP
0 sg
4030 3776 mt 3952 3195 L
3952 3195 mt 3868 3349 L
c26
-69 -608 -93 181 4030 3776 3 MP
PP
0 sg
4030 3776 mt 3937 3957 L
3937 3957 mt 3868 3349 L
/c27 { 0.812500 1.000000 0.187500 sr} bdef
c27
-88 -548 -74 160 4022 2310 3 MP
PP
0 sg
4022 2310 mt 3948 2470 L
3948 2470 mt 3860 1922 L
c27
-65 125 -97 -513 4022 2310 3 MP
PP
0 sg
4022 2310 mt 3925 1797 L
3925 1797 mt 3860 1922 L
c14
-45 65 -126 -386 3993 975 3 MP
PP
0 sg
3993  975 mt 3867  589 L
3867  589 mt 3822  654 L
c14
-116 -431 -55 110 3993 975 3 MP
PP
0 sg
3993  975 mt 3938 1085 L
3938 1085 mt 3822  654 L
c10
-98 -512 -64 136 3989 1661 3 MP
PP
0 sg
3989 1661 mt 3925 1797 L
3925 1797 mt 3827 1285 L
c10
-55 97 -107 -473 3989 1661 3 MP
PP
0 sg
3989 1661 mt 3882 1188 L
3882 1188 mt 3827 1285 L
/c28 { 0.000000 0.812500 1.000000 sr} bdef
c28
-68 -608 -84 154 3952 3195 3 MP
PP
0 sg
3952 3195 mt 3868 3349 L
3868 3349 mt 3800 2741 L
c28
-74 126 -78 -580 3952 3195 3 MP
PP
0 sg
3952 3195 mt 3874 2615 L
3874 2615 mt 3800 2741 L
/c29 { 0.562500 1.000000 0.437500 sr} bdef
c29
-79 -580 -74 145 3948 2470 3 MP
PP
0 sg
3948 2470 mt 3874 2615 L
3874 2615 mt 3795 2035 L
c29
-65 113 -88 -548 3948 2470 3 MP
PP
0 sg
3948 2470 mt 3860 1922 L
3860 1922 mt 3795 2035 L
c9
-46 61 -116 -431 3938 1085 3 MP
PP
0 sg
3938 1085 mt 3822  654 L
3822  654 mt 3776  715 L
c9
-106 -473 -56 103 3938 1085 3 MP
PP
0 sg
3938 1085 mt 3882 1188 L
3882 1188 mt 3776  715 L
/c30 { 0.000000 0.000000 1.000000 sr} bdef
c30
-83 125 -69 -608 3937 3957 3 MP
PP
0 sg
3937 3957 mt 3868 3349 L
3868 3349 mt 3785 3474 L
c30
-59 -630 -93 147 3937 3957 3 MP
PP
0 sg
3937 3957 mt 3844 4104 L
3844 4104 mt 3785 3474 L
/c31 { 1.000000 0.750000 0.000000 sr} bdef
c31
-88 -548 -65 125 3925 1797 3 MP
PP
0 sg
3925 1797 mt 3860 1922 L
3860 1922 mt 3772 1374 L
c31
-55 89 -98 -512 3925 1797 3 MP
PP
0 sg
3925 1797 mt 3827 1285 L
3827 1285 mt 3772 1374 L
/c32 { 1.000000 0.187500 0.000000 sr} bdef
c32
-97 -513 -55 97 3882 1188 3 MP
PP
0 sg
3882 1188 mt 3827 1285 L
3827 1285 mt 3730  772 L
c32
-46 57 -106 -473 3882 1188 3 MP
PP
0 sg
3882 1188 mt 3776  715 L
3776  715 mt 3730  772 L
/c33 { 0.312500 1.000000 0.687500 sr} bdef
c33
-64 99 -79 -580 3874 2615 3 MP
PP
0 sg
3874 2615 mt 3795 2035 L
3795 2035 mt 3731 2134 L
c33
-69 -607 -74 126 3874 2615 3 MP
PP
0 sg
3874 2615 mt 3800 2741 L
3800 2741 mt 3731 2134 L
/c34 { 0.000000 0.562500 1.000000 sr} bdef
c34
-75 104 -68 -608 3868 3349 3 MP
PP
0 sg
3868 3349 mt 3800 2741 L
3800 2741 mt 3725 2845 L
c34
-60 -629 -83 125 3868 3349 3 MP
PP
0 sg
3868 3349 mt 3785 3474 L
3785 3474 mt 3725 2845 L
/c35 { 1.000000 1.000000 0.000000 sr} bdef
c35
-55 81 -88 -548 3860 1922 3 MP
PP
0 sg
3860 1922 mt 3772 1374 L
3772 1374 mt 3717 1455 L
c35
-78 -580 -65 113 3860 1922 3 MP
PP
0 sg
3860 1922 mt 3795 2035 L
3795 2035 mt 3717 1455 L
/c36 { 0.000000 0.000000 0.750000 sr} bdef
c36
-49 -645 -94 106 3844 4104 3 MP
PP
0 sg
3844 4104 mt 3750 4210 L
3750 4210 mt 3701 3565 L
c36
-84 91 -59 -630 3844 4104 3 MP
PP
0 sg
3844 4104 mt 3785 3474 L
3785 3474 mt 3701 3565 L
c13
-46 54 -97 -513 3827 1285 3 MP
PP
0 sg
3827 1285 mt 3730  772 L
3730  772 mt 3684  826 L
c13
-88 -548 -55 89 3827 1285 3 MP
PP
0 sg
3827 1285 mt 3772 1374 L
3772 1374 mt 3684  826 L
/c37 { 0.125000 1.000000 0.875000 sr} bdef
c37
-59 -630 -75 104 3800 2741 3 MP
PP
0 sg
3800 2741 mt 3725 2845 L
3725 2845 mt 3666 2215 L
c37
-65 81 -69 -607 3800 2741 3 MP
PP
0 sg
3800 2741 mt 3731 2134 L
3731 2134 mt 3666 2215 L
c27
-69 -608 -64 99 3795 2035 3 MP
PP
0 sg
3795 2035 mt 3731 2134 L
3731 2134 mt 3662 1526 L
c27
-55 71 -78 -580 3795 2035 3 MP
PP
0 sg
3795 2035 mt 3717 1455 L
3717 1455 mt 3662 1526 L
/c38 { 0.000000 0.375000 1.000000 sr} bdef
c38
-74 75 -60 -629 3785 3474 3 MP
PP
0 sg
3785 3474 mt 3725 2845 L
3725 2845 mt 3651 2920 L
c38
-50 -645 -84 91 3785 3474 3 MP
PP
0 sg
3785 3474 mt 3701 3565 L
3701 3565 mt 3651 2920 L
/c39 { 1.000000 0.500000 0.000000 sr} bdef
c39
-45 49 -88 -548 3772 1374 3 MP
PP
0 sg
3772 1374 mt 3684  826 L
3684  826 mt 3639  875 L
c39
-78 -580 -55 81 3772 1374 3 MP
PP
0 sg
3772 1374 mt 3717 1455 L
3717 1455 mt 3639  875 L
/c40 { 0.000000 0.000000 0.562500 sr} bdef
c40
-84 47 -49 -645 3750 4210 3 MP
PP
0 sg
3750 4210 mt 3701 3565 L
3701 3565 mt 3617 3612 L
c40
-40 -651 -93 53 3750 4210 3 MP
PP
0 sg
3750 4210 mt 3657 4263 L
3657 4263 mt 3617 3612 L
/c41 { 0.687500 1.000000 0.312500 sr} bdef
c41
-59 -630 -65 81 3731 2134 3 MP
PP
0 sg
3731 2134 mt 3666 2215 L
3666 2215 mt 3607 1585 L
c41
-55 59 -69 -608 3731 2134 3 MP
PP
0 sg
3731 2134 mt 3662 1526 L
3662 1526 mt 3607 1585 L
/c42 { 0.000000 0.937500 1.000000 sr} bdef
c42
-50 -645 -74 75 3725 2845 3 MP
PP
0 sg
3725 2845 mt 3651 2920 L
3651 2920 mt 3601 2275 L
c42
-65 60 -59 -630 3725 2845 3 MP
PP
0 sg
3725 2845 mt 3666 2215 L
3666 2215 mt 3601 2275 L
/c43 { 1.000000 0.625000 0.000000 sr} bdef
c43
-69 -608 -55 71 3717 1455 3 MP
PP
0 sg
3717 1455 mt 3662 1526 L
3662 1526 mt 3593  918 L
c43
-46 43 -78 -580 3717 1455 3 MP
PP
0 sg
3717 1455 mt 3639  875 L
3639  875 mt 3593  918 L
/c44 { 0.000000 0.187500 1.000000 sr} bdef
c44
-74 40 -50 -645 3701 3565 3 MP
PP
0 sg
3701 3565 mt 3651 2920 L
3651 2920 mt 3577 2960 L
c44
-40 -652 -84 47 3701 3565 3 MP
PP
0 sg
3701 3565 mt 3617 3612 L
3617 3612 mt 3577 2960 L
c29
-56 45 -59 -630 3666 2215 3 MP
PP
0 sg
3666 2215 mt 3607 1585 L
3607 1585 mt 3551 1630 L
c29
-50 -645 -65 60 3666 2215 3 MP
PP
0 sg
3666 2215 mt 3601 2275 L
3601 2275 mt 3551 1630 L
c31
-46 38 -69 -608 3662 1526 3 MP
PP
0 sg
3662 1526 mt 3593  918 L
3593  918 mt 3547  956 L
c31
-60 -629 -55 59 3662 1526 3 MP
PP
0 sg
3662 1526 mt 3607 1585 L
3607 1585 mt 3547  956 L
c40
-84 -10 -40 -651 3657 4263 3 MP
PP
0 sg
3657 4263 mt 3617 3612 L
3617 3612 mt 3533 3602 L
c40
-30 -648 -94 -13 3657 4263 3 MP
PP
0 sg
3657 4263 mt 3563 4250 L
3563 4250 mt 3533 3602 L
/c45 { 0.000000 0.875000 1.000000 sr} bdef
c45
-65 33 -50 -645 3651 2920 3 MP
PP
0 sg
3651 2920 mt 3601 2275 L
3601 2275 mt 3536 2308 L
c45
-41 -652 -74 40 3651 2920 3 MP
PP
0 sg
3651 2920 mt 3577 2960 L
3577 2960 mt 3536 2308 L
c44
-75 -5 -40 -652 3617 3612 3 MP
PP
0 sg
3617 3612 mt 3577 2960 L
3577 2960 mt 3502 2955 L
c44
-31 -647 -84 -10 3617 3612 3 MP
PP
0 sg
3617 3612 mt 3533 3602 L
3533 3602 mt 3502 2955 L
c16
-45 29 -60 -629 3607 1585 3 MP
PP
0 sg
3607 1585 mt 3547  956 L
3547  956 mt 3502  985 L
c16
-49 -645 -56 45 3607 1585 3 MP
PP
0 sg
3607 1585 mt 3551 1630 L
3551 1630 mt 3502  985 L
/c46 { 0.500000 1.000000 0.500000 sr} bdef
c46
-40 -651 -65 33 3601 2275 3 MP
PP
0 sg
3601 2275 mt 3536 2308 L
3536 2308 mt 3496 1657 L
c46
-55 27 -50 -645 3601 2275 3 MP
PP
0 sg
3601 2275 mt 3551 1630 L
3551 1630 mt 3496 1657 L
c28
-31 -648 -75 -5 3577 2960 3 MP
PP
0 sg
3577 2960 mt 3502 2955 L
3502 2955 mt 3471 2307 L
c28
-65 -1 -41 -652 3577 2960 3 MP
PP
0 sg
3577 2960 mt 3536 2308 L
3536 2308 mt 3471 2307 L
/c47 { 0.000000 0.000000 0.625000 sr} bdef
c47
-84 -85 -30 -648 3563 4250 3 MP
PP
0 sg
3563 4250 mt 3533 3602 L
3533 3602 mt 3449 3517 L
c47
-21 -629 -93 -104 3563 4250 3 MP
PP
0 sg
3563 4250 mt 3470 4146 L
3470 4146 mt 3449 3517 L
c22
-40 -652 -55 27 3551 1630 3 MP
PP
0 sg
3551 1630 mt 3496 1657 L
3496 1657 mt 3456 1005 L
c22
-46 20 -49 -645 3551 1630 3 MP
PP
0 sg
3551 1630 mt 3502  985 L
3502  985 mt 3456 1005 L
c15
-55 3 -40 -651 3536 2308 3 MP
PP
0 sg
3536 2308 mt 3496 1657 L
3496 1657 mt 3441 1660 L
c15
-30 -647 -65 -1 3536 2308 3 MP
PP
0 sg
3536 2308 mt 3471 2307 L
3471 2307 mt 3441 1660 L
c26
-21 -628 -84 -85 3533 3602 3 MP
PP
0 sg
3533 3602 mt 3449 3517 L
3449 3517 mt 3428 2889 L
c26
-74 -66 -31 -647 3533 3602 3 MP
PP
0 sg
3533 3602 mt 3502 2955 L
3502 2955 mt 3428 2889 L
c45
-21 -628 -74 -66 3502 2955 3 MP
PP
0 sg
3502 2955 mt 3428 2889 L
3428 2889 mt 3407 2261 L
c45
-64 -46 -31 -648 3502 2955 3 MP
PP
0 sg
3502 2955 mt 3471 2307 L
3471 2307 mt 3407 2261 L
/c48 { 1.000000 0.937500 0.000000 sr} bdef
c48
-46 7 -40 -652 3496 1657 3 MP
PP
0 sg
3496 1657 mt 3456 1005 L
3456 1005 mt 3410 1012 L
c48
-31 -648 -55 3 3496 1657 3 MP
PP
0 sg
3496 1657 mt 3441 1660 L
3441 1660 mt 3410 1012 L
c46
-21 -629 -64 -46 3471 2307 3 MP
PP
0 sg
3471 2307 mt 3407 2261 L
3407 2261 mt 3386 1632 L
c46
-55 -28 -30 -647 3471 2307 3 MP
PP
0 sg
3471 2307 mt 3441 1660 L
3441 1660 mt 3386 1632 L
/c49 { 0.000000 0.000000 0.875000 sr} bdef
c49
-12 -588 -92 -236 3470 4146 3 MP
PP
0 sg
3470 4146 mt 3378 3910 L
3378 3910 mt 3366 3322 L
c49
-83 -195 -21 -629 3470 4146 3 MP
PP
0 sg
3470 4146 mt 3449 3517 L
3449 3517 mt 3366 3322 L
/c50 { 0.000000 0.437500 1.000000 sr} bdef
c50
-12 -587 -83 -195 3449 3517 3 MP
PP
0 sg
3449 3517 mt 3366 3322 L
3366 3322 mt 3354 2735 L
c50
-74 -154 -21 -628 3449 3517 3 MP
PP
0 sg
3449 3517 mt 3428 2889 L
3428 2889 mt 3354 2735 L
c22
-46 -8 -31 -648 3441 1660 3 MP
PP
0 sg
3441 1660 mt 3410 1012 L
3410 1012 mt 3364 1004 L
c22
-22 -628 -55 -28 3441 1660 3 MP
PP
0 sg
3441 1660 mt 3386 1632 L
3386 1632 mt 3364 1004 L
c18
-12 -588 -74 -154 3428 2889 3 MP
PP
0 sg
3428 2889 mt 3354 2735 L
3354 2735 mt 3342 2147 L
c18
-65 -114 -21 -628 3428 2889 3 MP
PP
0 sg
3428 2889 mt 3407 2261 L
3407 2261 mt 3342 2147 L
/c51 { 0.625000 1.000000 0.375000 sr} bdef
c51
-55 -73 -21 -629 3407 2261 3 MP
PP
0 sg
3407 2261 mt 3386 1632 L
3386 1632 mt 3331 1559 L
c51
-11 -588 -65 -114 3407 2261 3 MP
PP
0 sg
3407 2261 mt 3342 2147 L
3342 2147 mt 3331 1559 L
c16
-45 -32 -22 -628 3386 1632 3 MP
PP
0 sg
3386 1632 mt 3364 1004 L
3364 1004 mt 3319  972 L
c16
-12 -587 -55 -73 3386 1632 3 MP
PP
0 sg
3386 1632 mt 3331 1559 L
3331 1559 mt 3319  972 L
c38
-81 -389 -12 -588 3378 3910 3 MP
PP
0 sg
3378 3910 mt 3366 3322 L
3366 3322 mt 3285 2933 L
c38
-3 -508 -90 -469 3378 3910 3 MP
PP
0 sg
3378 3910 mt 3288 3441 L
3288 3441 mt 3285 2933 L
c45
-3 -509 -81 -389 3366 3322 3 MP
PP
0 sg
3366 3322 mt 3285 2933 L
3285 2933 mt 3282 2424 L
c45
-72 -311 -12 -587 3366 3322 3 MP
PP
0 sg
3366 3322 mt 3354 2735 L
3354 2735 mt 3282 2424 L
/c52 { 0.375000 1.000000 0.625000 sr} bdef
c52
-63 -231 -12 -588 3354 2735 3 MP
PP
0 sg
3354 2735 mt 3342 2147 L
3342 2147 mt 3279 1916 L
c52
-3 -508 -72 -311 3354 2735 3 MP
PP
0 sg
3354 2735 mt 3282 2424 L
3282 2424 mt 3279 1916 L
/c53 { 0.875000 1.000000 0.125000 sr} bdef
c53
-55 -152 -11 -588 3342 2147 3 MP
PP
0 sg
3342 2147 mt 3331 1559 L
3331 1559 mt 3276 1407 L
c53
-3 -509 -63 -231 3342 2147 3 MP
PP
0 sg
3342 2147 mt 3279 1916 L
3279 1916 mt 3276 1407 L
c43
-3 -508 -55 -152 3331 1559 3 MP
PP
0 sg
3331 1559 mt 3276 1407 L
3276 1407 mt 3273  899 L
c43
-46 -73 -12 -587 3331 1559 3 MP
PP
0 sg
3331 1559 mt 3319  972 L
3319  972 mt 3273  899 L
c33
-60 -981 -3 -508 3288 3441 3 MP
PP
0 sg
3288 3441 mt 3285 2933 L
3285 2933 mt 3225 1952 L
c33
1 -327 -64 -1162 3288 3441 3 MP
PP
0 sg
3288 3441 mt 3224 2279 L
3224 2279 mt 3225 1952 L
c51
-55 -799 -3 -509 3285 2933 3 MP
PP
0 sg
3285 2933 mt 3282 2424 L
3282 2424 mt 3227 1625 L
c51
2 -327 -60 -981 3285 2933 3 MP
PP
0 sg
3285 2933 mt 3225 1952 L
3225 1952 mt 3227 1625 L
c35
1 -326 -55 -799 3282 2424 3 MP
PP
0 sg
3282 2424 mt 3227 1625 L
3227 1625 mt 3228 1299 L
c35
-51 -617 -3 -508 3282 2424 3 MP
PP
0 sg
3282 2424 mt 3279 1916 L
3279 1916 mt 3228 1299 L
/c54 { 1.000000 0.687500 0.000000 sr} bdef
c54
-46 -435 -3 -509 3279 1916 3 MP
PP
0 sg
3279 1916 mt 3276 1407 L
3276 1407 mt 3230  972 L
c54
2 -327 -51 -617 3279 1916 3 MP
PP
0 sg
3279 1916 mt 3228 1299 L
3228 1299 mt 3230  972 L
c23
-42 -254 -3 -508 3276 1407 3 MP
PP
0 sg
3276 1407 mt 3273  899 L
3273  899 mt 3231  645 L
c23
1 -327 -46 -435 3276 1407 3 MP
PP
0 sg
3276 1407 mt 3230  972 L
3230  972 mt 3231  645 L
gr


end %%Color Dict

eplot
%%EndObject

396 18 moveto
/Helvetica-BoldOblique findfont
[10 0 0 10 0 0] makefont setfont
(Student Version of MATLAB) show

epage
end

showpage

%%Trailer
%%BoundingBox:    66   212   567   583
%%Pages: 001
%%EOF
