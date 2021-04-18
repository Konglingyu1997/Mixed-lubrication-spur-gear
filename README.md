From file: Predicting_Spur gear0418_Only simulate single contact unit
正齒輪
設計參數

m = 4;(*模數,mm*)
\[Alpha] = 20 Degree;(*法面壓力角*)
z1 = 21;(*主動輪G1的齒數*)
z2 = 49;(*從動輪G2的齒數*)
\[Omega]1 = N[ 2, 16 ];(*轉速，rad/s*)
B = 20;(*face width,mm*)
dn = 2.5*10^-3;(*粗糙面節點間距,mm*)


二級參數

齒輪幾何參數

(*數字1、2分別表示gear1和gear2*)
\[Omega]2 = N[z1/z2*\[Omega]1, 16];
rp1 = N[z1*m/2, 16];
rp2 = N[z2*m/2, 16];
(*基圓柱半徑*)
rb1 = N[rp1*Cos[\[Alpha]], 16];
rb2 = N[rp2*Cos[\[Alpha]], 16];
(*齒頂圓柱半徑*)
ra1 = N[rp1 + m, 16];
ra2 = N[rp2 + m, 16];
(*端面基節,base pitch*)
Pb = N[2*Pi*rb1/z1, 16];
(*接觸率*)
(*length of contact line at transverse plane*)
Lc = N[Sqrt[ra1^2 - rb1^2] + Sqrt[ra2^2 - rb2^2] - (rp1 + rp2)*Sin[\[Alpha]], 
   16];
\[Epsilon]\[Alpha] = N[ Lc/Pb, 16 ];(*transverse contact ratio*)
(*Fig.1(b)各條線段長度*)
N1N2 = N[(rb1 + rb2)*Tan[\[Alpha]], 16];
N2A = N[Sqrt[ra2^2 - rb2^2], 16];
N1A = N[N1N2 - N2A, 16];(*式[1f]Subscript[中之lN, 1]A*)

動力學參數

Fn = 500;(*transmitted force,N*)
eqE = 2.07*10^11;(*equivalent elastic module,N/m^2*)
eqE = eqE*10^-6;(*N/mm^2*)


流體力學參數

\[Eta]0 = 0.08;(*  Pa*s  *)
\[Eta]0 = \[Eta]0*10^-6;(*   N/mm^2*s   *)


求解器設定

(*求解器設定*)
(*總時間步數*)
stepnum = 30;
(*作用線上節點數量*)
ynum = 30;
(*有效接觸齒寬上節點數量*)
unum = 20;
(*三個特殊時間節點*)
(*t1=N[ B*Tan[\[Beta]b]/(rb1*\[Omega]1),16 ];
t2=N[ Lc/(rb1*\[Omega]1),16 ];
t3=N[ Pbt1*(\[Epsilon]\[Alpha]+\[Epsilon]\[Beta])/(rb1*\[Omega]1),16 ];*)
(*單齒嚙合時間*)
t3 = N[Lc/(rb1*\[Omega]1), 16]
(*單位時間步長*)
dt = N[t3/stepnum, 16]
(*在接觸線上的單位距離間隔,單位：mm*)
dl = N[B/unum, 16];
(*在有效接觸齒寬（接觸線在齒寬方向之投影）上的單位距離間隔,單位：mm*)
du = N[B/unum, 16];


0.248356188470812

0.0082785396156937


Dynamic analysis
(*u為有效嚙合齒寬上的位置參數*)
Rxstage1 = {};(*equuivalent radius of curvature*)
Uxstage1 = {};(*rolling velocity*)
Uxsstage1 = {};(*relative sliding velocity*)
For[t = 0, t <= t3, t += dt,
  
  For[u = 0, u <= B, u += du,
    (*等效曲率半徑*)
    R1x[t][u] = N1A + rb1*\[Omega]1*t + u*0;
    R2x[t][u] = N1N2 - R1x[t][u];
    Rx[t][u] = ( R1x[t][u]*R2x[t][u]/(R1x[t][u] + R2x[t][u]) );
    AppendTo[  Rxstage1, {t, u, Rx[t][u]} ];
    (*速度*)
    U1x[t][u] = \[Omega]1*R1x[t][u];
    U2x[t][u] = \[Omega]2*R2x[t][u];
    (*rolling velocity,即entrainment velocity*)
    Ux[t][u] = (U1x[t][u] + U2x[t][u])/2;
    AppendTo[  Uxstage1, {t, u, Ux[t][u]} ];
    (*relative sliding velocity*)
    Uxs[t][u] = (U1x[t][u] - U2x[t][u]);
    AppendTo[  Uxsstage1, {t, u, Uxs[t][u]} ];
    ];
  
  ];
  
  
  
  
  (*嚙入，全齒寬嚙合，嚙出，三個階段等效曲率半徑 圖像。
x軸：時間t，單位s；
y軸：接觸線上之位置參數u，單位：mm；
z軸：等效曲率半徑Rx，單位：mm*)
size = 12;
labelx = Text[ Style[ "meshing time(s)", FontSize -> size] ];
labely = Text[ Style["tooth width(mm)", FontSize -> size] ];
labelz = Text[ Style["Rx(mm)", FontSize -> size] ];
pRxstage1 = 
 ListPointPlot3D[ Rxstage1, PlotStyle -> Orange, 
  AxesLabel -> {labelx, labely, labelz}, 
  LabelStyle -> Directive[Black, Small, Bold], PlotRange -> All, 
  Boxed -> False ];
  
 (*嚙入，全齒寬嚙合，嚙出，三個階段相對滑動速度 圖像。
x軸：時間t，單位s；
y軸：接觸線上之位置參數u，單位：mm；
z軸：相對滑動速度Uxs，單位：mm/s*)
size = 12;
labelx = Text[ Style[ "meshing time(s)", FontSize -> size] ];
labely = Text[ Style["tooth width(mm)", FontSize -> size] ];
labelz = Text[ Style["Uxs(mm/s)", FontSize -> size] ];
pUxsstage1 = 
  ListPointPlot3D[ Uxsstage1, PlotStyle -> Orange, 
   AxesLabel -> {labelx, labely, labelz}, 
   LabelStyle -> Directive[Black, Small, Bold] ];
Show[pUxsstage1, PlotRange -> All, Boxed -> False ]


Mechanical analysis
Load sharing ratio

\[Xi]inn = N[N1A/Pb, 16](*profile parameter at inner point*)
\[Xi]a = N[(Lc + N1A)/Pb, 16](*profile parameter at aduddenm point*)

t = 0;
While[True,
  \[Xi][t] = (\[Omega]1*t*rb1 + N1A)/Pb;
  If[ \[Xi][t] > \[Xi]a - 1,
   \[Xi][t] = Null;
   Break[] 
   ];
  R[t] = 1/3*(  1 + (\[Xi][t] - \[Xi]inn)/(\[Xi]a - \[Xi]inn - 1)  );
  t += dt;
  ];
t1 = t - dt;

While[True,
  \[Xi][t] = (\[Omega]1*t*rb1 + N1A)/Pb;
  If[ \[Xi][t] > \[Xi]inn + 1,
   \[Xi][t] = Null;
   Break[] 
   ];
  R[t] = 1;
  t += dt;
  ];
t2 = t - dt;

While[True,
  \[Xi][t] = (\[Omega]1*t*rb1 + N1A)/Pb;
  If[ \[Xi][t] > \[Xi]a,
   \[Xi][t] = Null;
   Break[] 
   ];
  R[t] = 1/3*(  1 + (\[Xi][t] - \[Xi]a)/(\[Xi]inn + 1 - \[Xi]a)  );
  t += dt;
  ];
t3 = t - dt;
(*load sharing ratio(R) for spur gears*)
data = {};
For[t = 0, t <= t3, t += dt,
  AppendTo[data, {t, R[t]}];
  ];
labelx = Text[  Style[ "t(s)", FontSize -> 20 ]  ];
labely = Text[  Style["R", FontSize -> 20]  ];
fig3 = ListPlot[ data, AxesLabel -> { labelx , labely  }, 
  LabelStyle -> Directive[Black, Bold] ]
labelx = Text[ Style["\[Xi]", FontSize -> 20] ];
labelx1 = 
  Text[ Style["profile parameter on transverse plane", 
    FontSize -> 14] ];
labely = Text[ Style["R", FontSize -> 20] ];
labely1 = Text[ Style["load sharing ratio", FontSize -> 14] ];
fig4 = ListPlot[  Table[{\[Xi][t], R[t]}, {t, 0, t3, dt}], 
  AxesLabel -> {labelx, labely}, LabelStyle -> Directive[Black, Bold],
   Frame -> True, FrameLabel -> {labelx1, labely1}];

Load distribution

data1 = {};
For[t = 0, t <= t3, t += dt,
  For[u = 0, u <= B, u += du,
    f[t][u] = Fn*R[t]/B + u*0;
    AppendTo[  data1, {t, u, f[t][u]} ];
    ];
  ];
size = 16;(*字體大小*)
labelx = Text[ Style[ "meshing time(s)", FontSize -> size] ];
labely = Text[ 
   Style["length on line of action(mm)", FontSize -> size] ];
labelz = Text[ Style["f(N/mm)", FontSize -> size] ];
fig4 = ListPointPlot3D[ data1, 
  AxesLabel -> { labelx, labely, labelz }, 
  LabelStyle -> Directive[Black, Medium, Bold] ];
 
 
 Contact analysis
 
 
 (* Maximum Hertzian pressure ph at each node
定義
x方向：漸開線齒廓方向（profile direction)
y方向：齒長方向（axial direction） *)

(*赫茲接觸中心：圖中負載W正下方的點，即瞬時接觸坐標之y軸經過的點
以各節點作為瞬時的赫茲接觸中心，計算各節點的最大赫茲壓力。
力F: 單位長度之受力乘兩點間距；
接觸線長度B: 在齒長方向(y方向)之相鄰兩節點間距du；
半徑Rx: x方向的等效曲率半徑*)

data2 = {};
data3 = {};(*Note!!!data3/data4中儲存的是赫茲接觸寬度，不是半寬*)
data4 = {};
For[t = dt, t <= t3 - dt, t += dt,(*嚙入、嚙出點為半球接觸，赫茲壓力無法計算,故距離起始於dt,
  終了於t3-dt*)
  For[u = 0, u <= B, u += du,
    (*各點所受正向力(沿接觸線方向之受力);單位：F,N;  f,N/mm;  du,mm*)
    F[t][u] = f[t][u]*du;
    (*各點最大赫茲壓力;單位：ph,N;  F,N;  eqE,N/mm^2;  du,mm;  Rx,mm*)
    ph[t][u] = ( F[t][u]*eqE/(2*Pi*du*Rx[t][u])  )^(1/2);
    (*赫茲接觸半寬;單位：a,mm;  F,N;  Rx,mm;  eqE,N/mm^2*)
    a[t][u] = ( 3*F[t][u]*Rx[t][u]/(4*eqE) )^(1/3);
    AppendTo[ data2, {t, u, ph[t][u]} ];
    (*Note!!!data3/data4中儲存的是赫茲接觸寬度，不是半寬*)
    AppendTo[ data3, {t, u, 2*a[t][u]} ];
    AppendTo[ data4, 2*a[t][u] ];
    ];
  ];
size = 12;(*字體大小*)
labelx = Text[ Style[ "meshing time(s)", FontSize -> size] ];
labely = Text[ 
   Style["length on line of action(mm)", FontSize -> size] ];
labelz = Text[ 
   Style["maximum Hertz pressure(N/\!\(\*SuperscriptBox[\(mm\), \
\(2\)]\))", FontSize -> size] ];
fig4 = ListPointPlot3D[ data2, 
  AxesLabel -> { labelx, labely, labelz }, 
  LabelStyle -> Directive[Black, Tiny, Bold] , Boxed -> False];
  
  
 size = 12;(*字體大小*)
labelx = Text[ Style[ "meshing time(s)", FontSize -> size] ];
labely = Text[ Style["length(mm)", FontSize -> size] ];
labelz = Text[ Style["HHCW(mm)", FontSize -> size] ];
fig4 = ListPointPlot3D[ data3, 
  AxesLabel -> { labelx, labely, labelz }, 
  LabelStyle -> Directive[Black, Small, Bold] , Boxed -> False];
  
 
 
 (*Hertz contact interval at each node
定義：1.作用線（line of action）：x軸
2.接觸線（line of contact）：y軸
赫茲接觸區間：設接觸中心在plane of action上的坐標為(x1,y1)，接觸半寬為a，則赫茲接觸區間為[x1-a,x1+a]
定理：從嚙入點到嚙合過程中之某時刻，作用點在齒廓上經過的漸開線弧長，等於瞬時作用線段長
步驟1：求得各節點在plane of action \
上的橫坐標x1，即各節點之瞬時坐標係（x'Oy'）的原點O，在全域坐標x軸上的位置，記為con[t]
步驟2：求得各節點對應之赫茲接觸區間，區間左端點為xle[t]=x1-a，右端點為xri[t]=x1+a;*)
(*步驟1、2*)
data5 = {};
data6 = {};
data7 = {};
For[t = dt, t <= t3 - dt, t += dt,
  con[t] = rb1*\[Omega]1*t;
  xle[t] = rb1*\[Omega]1*t - a[t][0];
  xri[t] = rb1*\[Omega]1*t + a[t][0];
  AppendTo[ data5, con[t] ];
  AppendTo[ data6, xle[t] ];
  AppendTo[ data7, xri[t] ];
  ];
  
  (*步驟3：構造無因次接觸坐標軸，並求左右端點對應的首週期無因次坐標
A.定義：
B. 定義x'軸截得之表面高度曲線
首先，觀察plane of action，在接觸線上，各接觸點之縱坐標沿接觸線方向等距分佈，距離du=1mm
又因為每個磨紋單元之邊長為0.25mm，所以相鄰兩個接觸點在y方向之距離=磨紋單元邊長*4。
而接觸點O在y方向之坐標，決定了每個瞬時接觸坐標系之x'軸所截得之表面高度曲線，所以對任意接觸點，其x'\
軸截得之表面高度曲線相同。此案例定義為磨紋單元的前端截面之表面曲線
C. 定義無因次接觸坐標軸：
方向：以作用線方向為x軸，作用點移動方向為正方向
節點間距dn=2.5*10^-3mm，單個週期節點數量T=100，單週期長度LT為單個磨紋單元之邊長，LT=T*dn
節點編號為i，最左端i=0，最右端i=99; 
記接觸區左端點在磨紋單元中對應的節點為nle[t]；右端點對應nri[t]
取點原則：當接觸區端點與節點不重合時，左端點向右取，右端點向左取，即盡量縮小接觸區求
1）各區間之左右端點，分別在粗糙度函數之無因次坐標軸上對應的點，
2）並利用其週期性，求函數值;
3）若所有點中最大接觸寬度2a大於週期T之長度，求接觸區包含的函數週期數（此案不用）*)

dn = 2.5*10^-3;(*mm*)
T = 100;(*x方向的粗糙度函數之週期*)
For[ t = dt, t <= t3 - dt, t += dt,
  (*t時刻，接觸區左端點對應的無因次坐標*)
  x1[t] = Ceiling[   xle[t]/dn  ];
  (*接觸區左端點對應的首週期節點坐標*)
  n1[t] = Mod[  x1[t]  , T  ];
  (*t時刻，接觸區右端點對應的無因次坐標*)
  x2[t] = Floor[  xri[t]/dn  ];
  (*接觸區右端點對應的首週期節點坐標*)
  n2[t] = Mod[  x2[t]  , T  ];
  ];

(*各接觸區：1. 無因次寬度
2. 實際寬度
3. 無因次寬度之誤差(無因次寬度-實際寬度)*)

dimenwid = ListPlot[  Table[x2[t] - x1[t], {t, dt, t3 - dt, dt}]  ];
realwid = ListPlot[  Table[2*a[t][0], {t, dt, t3 - dt, dt}]  ];
error = ListPlot[ 
   Table[ x2[t] - x1[t] - 2*a[t][0]/dn, {t, dt, t3 - dt, dt}] ];

  
