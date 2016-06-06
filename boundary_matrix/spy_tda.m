function plot_handle = spy_tda(S,arg2,arg3)
%SPY Visualize sparsity pattern.
%   SPY(S) plots the sparsity pattern of the matrix S.
%
%   SPY(S,'LineSpec') uses the color and marker from the line
%   specification string 'LineSpec' (See PLOT for possibilities).
%
%   SPY(S,markersize) uses the specified marker size instead of
%   a size which depends upon the figure size and the matrix order.
%
%   SPY(S,'LineSpec',markersize) sets both.
%
%   SPY(S,markersize,'LineSpec') also works.

%   Copyright 1984-2013 The MathWorks, Inc. 
%   Modified May 2016 to return handles

cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

marker = ''; color = ''; markersize = 0; linestyle = 'none';
if nargin >= 2
   if ischar(arg2), 
      [~,color,marker,msg] = colstyle(arg2);
      if ~isempty(msg)
         error(message(msg.identifier));
      end
   else
      markersize = arg2;
   end
end
if nargin >= 3
   if ischar(arg3),
      [~,color,marker,msg] = colstyle(arg3); 
      if ~isempty(msg)
         error(message(msg.identifier));
      end
   else
      markersize = arg3;
   end
end
if isempty(marker), marker = '.'; end
if isempty(color), co = get(cax,'colororder'); color = co(1,:); end

if nargin < 1, S = defaultspy; end
[m,n] = size(S);
if marker~='.' && markersize==0,
   markersize = get(gcf,'defaultlinemarkersize');
end
if markersize == 0
   units = get(gca,'units');
   set(gca,'units','points');
   pos = get(gca,'position');
   markersize = max(4,min(14,round(6*min(pos(3:4))/max(m+1,n+1))));
   set(gca,'units',units);
end

[i,j] = find(S);
if isempty(i), i = NaN; j = NaN; end
if isempty(S), marker = 'none'; end
plot_handle = plot(j,i,'marker',marker,'markersize',markersize, ...
   'linestyle',linestyle,'color',color);

xlabel(['nz = ' int2str(nnz(S))]);
set(gca,'xlim',[0 n+1],'ylim',[0 m+1],'ydir','reverse', ...
   'GridLineStyle','none','plotboxaspectratio',[n+1 m+1 1]);

if ~hold_state, set(cax,'NextPlot',next); end

% ------------------------------------------------------

function S = defaultspy
c = [
    'a:321mjYCCXYFH9oo29:;bcdd;f7lz{`[!C_zYQNVzLORVyNVVV_]xL,1y5{92:y7>?>cL@['
    'b;796nnZYYYmGI:qr3:<?heeedg8m''!a\#Ea{{ROW{MPSWzOWWW`^yM-2z6*:3;z8?@?fMAf'
    'c_8::.o[ZZmqHJ;r*4;@@ifffej9r("b^=Fb?GSPXINQTX{PXxXwvzW.3{7+;4<{9@A@gNBy'
    'dd<<;/-kllorIK<*+9<A`jgggfkh''t#y_>Gc@ITQYJORUYJQyyYxw(X/4)80i5=3:ABAxOC:'
    '_f===0.mnmp*RLA+0:A^akhhhgli(x${`@HzA{URZKPSVZKRzzxyx)Y05/91j6f4;BCByPDF'
    '`g]>>5/oooq+SME03;\_hliiijpumy%?aAY{B{VSzLQTWbLS{{yzy*w1m0k2z7g5<ODCzQEG'
    'ah^\?;3ppp+0TNF34A]`imjjkku(uzx@dEZ>C{zT{MRUXyMTMQz{z+x2w1l3{8y6=PED{RFH'
    'b4c][>4-,q,?UOG49Z^hjokkllv3(&yA{F[?D{{UINSVYzNUNR{F{Iy3x2x409z7>QFE9SGI'
    'c5df\?<.-+0CVPH9:[iikpllpp(4,''{B!K\@E{MVJOTWZ{OVOSDG)Jz4y3y51:{8?SGF:THJ'
    'd;fj]@@21,>DWQI:;\jjlqppqt25nw!C>Z]AF{NWKPUX[JPWPTELGK{5z4z62;29@TIG;UIK'
    'e<ikg[A=21@EXRJ;Adkkosqqsu36uz"]?[cBG{OYLQVYbKQXQUFMHL''O{5{73<3:AbOH<VJL'
    '5^jlh\ZA=>CWYSKAYhmlptsst''4iv%#^@\dCH{PZMRWZyLRYRVMNIR,P(6*84=4;BcPM=WKM'
    '6_20ii[B>BDXmTLIdioos)ttu(5j(&$yA]zDI{QzNSX[zMSySWNONS-Q)7/95g5<cfQN>XLN'
    '7g31lmjYBCXYqUMJhjpp)*))()6r)wy{B_{EY{R{OTYb{NTzTZOPOTMR-80:6h6=dgSO?YMO'
    '8h8:mnkZYYYmrVNKikrr*1**)*7u3zz"C`=F[{SJPUZyJOU{U[PQPZNS.l1j7y7>fxTP@ZNP'
    '939</.mkZZmp*WOLjmss1211*19(4#z?Da>G^{TKQV[zKPV?V\QRQ[OT/m2k8z8?gy`QA[OQ'
    ':4=>0/nllmqq+XPMkr))233312:65$!BEc?H_yULRWb{LQWLWyRSRwSU0y3{9{9dy3aRB\PS'
    '_<>?64onoo*r,mQNms**344423g7w%zDFd@Iz{VMSXw<MRXMyzSTSxTl1z40:1:e44bSCcQT'
    '`=]\;5-opp+*0nROn)01455534h8)w!\[zEY{{WNTYzJNSYNz{TUTyUv2{51;2;g67cTDfRU'
    'd]^]=;.pqq0+:oSPo*22566645uj2xz_]{F\A{XOUZ{KOTZO{KUVUzVw3)62h3<y78gUEySV'
    'e^ef?<3,++>0<qTQr+336:8756(k3z!aa=G]B{YPV[ILPUyPPLVWV''vx4/73i4=089xWF:TW'
    'f_fg[><-,,?<@rVis044:;:8672qu{ycc>H^C{zQWzJMQVzQQQWXW(wy5084y5>19:9XG;UX'
    'gchh\@=2-1@@D)Wj)255;<;9783tw"zdd?I_D{{RX{KNRW{RRRXYX+xz6195z6?3:;:YH<WY'
    '4hjkfZ@=1>CCE*ik*366>=<:894ux#!{e@\`E{KSYILOSXKSSSYZ[,y(72k6{7@4;<;]I=XZ'
    '5iklh[AA>CDDF+jm+477?>=;9:5(y$^=zA]aF{LTZJMPTYLTTTZ[\Jz/m3y708e5<=<^J>Y['
    '6210iiBBBDXEG0nn069:abc<:;6k(x_D{B^bIPMU[KNQUZMUUU^\wK''0x4z819f6=>=_K?Z]'
    ];
i = double(c(:)-32);
j = cumsum(diff([0; i])<=0) + 1;
S = sparse(i,j,1)';
