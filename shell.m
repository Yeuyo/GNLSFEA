function [ke,fint,epsE,Vn,sig,L]=shell(nodeData,uvw,duvw,epsEn,oVn,...
                                       E,nu,ka,sigN,oL)
coord=nodeData(:,1:3); t=nodeData(:,6); epsE=zeros(6,8);
ke=zeros(9*5); fint=zeros(9*5,1); dxr=zeros(3); [wp,GpLoc]=GpPos();
ex=[1 0 0].'; ey=[0 1 0].'; ez=[0 0 1].';
V1=zeros(9,3); V2=V1; g1=V1; g2=V1;
xsi=[-1; -1; -1; 0; 1; 1; 1; 0; 0]; eta=[-1; 0; 1; 1; 1; 0; -1; -1; 0];
dNr=dershapefunc(xsi,eta); dnr=dNr(1:2:end,:); dns=dNr(2:2:end,:);
dsp=reshape(uvw-duvw,5,[])'; elcoord=coord+dsp(:,1:3);
for n=1:9
  Vn(n,:)=cross((dnr(n,:)*elcoord)/norm(dnr(n,:)*elcoord),...
                (dns(n,:)*elcoord)/norm(dns(n,:)*elcoord));
  V=cross(ey,Vn(n,:).').'; 
  V1(n,:)=V/norm(V); V2(n,:)=cross(Vn(n,:).',V1(n,:).').';
  g1(n,:)=-0.5*t(n)*V2(n,:); g2(n,:)= 0.5*t(n)*V1(n,:);
end
D=[[1 nu; nu 1; 0 0] zeros(3,4); zeros(3) eye(3)*ka*(1-nu)*0.5];
D=(E/(1-nu^2))*D; D(4,4)=D(4,4)/ka;
for Gp=1:8
  xsi=GpLoc(Gp,1); eta=GpLoc(Gp,2); zet=GpLoc(Gp,3);
  N=shapefunc(xsi,eta); dNr=dershapefunc(xsi,eta);
  dxr(1:2,:)=(dNr*(coord+0.5*zet*(oVn.*(t*ones(1,3)))));
  dxr(3,:)=(0.5*(N.*t')*oVn); VnGp=(N*oVn).';
  sdir=dxr(:,2)/norm(dxr(:,2));
  er=cross(sdir,VnGp); er=er/norm(er);
  es=cross(VnGp,er);   es=es/norm(es);
  et=VnGp/norm(VnGp);
  l1=(ex.'*er); m1=(ey.'*er); n1=(ez.'*er);
  l2=(ex.'*es); m2=(ey.'*es); n2=(ez.'*es);
  l3=(ex.'*et); m3=(ey.'*et); n3=(ez.'*et);
  Q=[l1*l1     m1*m1   n1*n1 l1*m1       m1*n1       n1*l1       ;
     l2*l2     m2*m2   n2*n2 l2*m2       m2*n2       n2*l2       ; 
     l3*l3     m3*m3   n3*n3 l3*m3       m3*n3       n3*l3       ; 
     2*l1*l2 2*m1*m2 2*n1*n2 l1*m2+l2*m1 m1*n2+m2*n1 n1*l2+n2*l1 ;
     2*l2*l3 2*m2*m3 2*n2*n3 l2*m3+l3*m2 m2*n3+m3*n2 n2*l3+n3*l2 ;
     2*l3*l1 2*m3*m1 2*n3*n1 l3*m1+l1*m3 m3*n1+m1*n3 n3*l1+n1*l3];
  Dsh=Q.'*D*Q;
  [B,GNL,epsEt,L(Gp,:),detJ]=formB(GpLoc(Gp,:),coord,oVn,t,g1,g2,...
                                   duvw,oL(Gp,:));
  epsE(:,Gp)=epsEn(:,Gp)+epsEt; sig(:,Gp)=sigN(:,Gp)+(Dsh*epsEt);
  H=[sig(1,Gp)*eye(3) sig(4,Gp)*eye(3) sig(6,Gp)*eye(3);
     sig(4,Gp)*eye(3) sig(2,Gp)*eye(3) sig(5,Gp)*eye(3);
     sig(6,Gp)*eye(3) sig(5,Gp)*eye(3) sig(3,Gp)*eye(3)];
  ke=ke+((B.'*Dsh*B)+(GNL.'*H*GNL))*detJ*wp(Gp);
  fint=fint+B.'*sig(:,Gp)*detJ*wp(Gp);
end

function [wp,GpLoc]=GpPos()
wp=ones(8,1); g2=1/sqrt(3); xsi=[-1 -1 1 1 -1 -1 1 1].'*g2;
eta=[-1 -1 -1 -1 1 1 1 1].'*g2; zet=[-1 1 1 -1 -1 1 1 -1].'*g2;
GpLoc=[xsi eta zet];

function [B,GNL,epsEt,L,detJ]=formB(GpLoc,coord,Vn,t,g1,g2,duvw,oL)
xsi=GpLoc(1); eta=GpLoc(2); zet=GpLoc(3);
N=shapefunc(xsi,eta); dNr=dershapefunc(xsi,eta); dxr=zeros(3);
dxr(1:2,:)=(dNr*(coord+0.5*zet*(Vn.*(t*ones(1,3)))));
dxr(3,:)=(0.5*(N.*t')*Vn);
detJ=det(dxr); invJ=inv(dxr); dNx=dxr\[dNr; zeros(1,9)];
G=zet*(invJ(:,1:2)*dNr)+invJ(:,3)*N;
B9=zeros(9,9*5); BNL=zeros(6,9*5); GNL=B9; L=zeros(1,9);
B9([1 4 9],1:5:end)=dNx; B9([5 2 6],2:5:end)=dNx;
B9([8 7 3],3:5:end)=dNx;
B9(:,4:5:end)=g1(:,[1 2 3 1 2 2 3 3 1])'.*G([1 2 3 2 1 3 2 1 3],:);
B9(:,5:5:end)=g2(:,[1 2 3 1 2 2 3 3 1])'.*G([1 2 3 2 1 3 2 1 3],:);
BL=B9([1:3 5 7 9],:); BL(4:6,:)=BL(4:6,:)+B9([4 6 8],:);
for n=1:3
  ct=3*(n-1)+1:3*(n-1)+3;
  for m=1:9
    m5=5*(m-1);
    L(:,ct)=L(:,ct)+([dNx(:,m) g1(m,n).*ones(3,1).*G(:,m)...
                      g2(m,n).*ones(3,1).*G(:,m)]*[duvw((m5)+n);
                      duvw((m5)+4); duvw((m5)+5)])';
  end
end
tL=oL+L;
A=[tL(:,1:3:9)  zeros(1,3)  zeros(1,3);
    zeros(1,3) tL(:,2:3:9)  zeros(1,3);
    zeros(1,3)  zeros(1,3) tL(:,3:3:9);
   tL(:,2:3:9) tL(:,1:3:9)  zeros(1,3);
    zeros(1,3) tL(:,3:3:9) tL(:,2:3:9);
   tL(:,3:3:9)  zeros(1,3) tL(:,1:3:9)];
for n=1:9
  ct=5*(n-1)+1:5*(n-1)+5;
  GNL(:,ct)=[dNx(1,n)*eye(3) g1(n,:)'.*G(1,n).*ones(3,1)...
             g2(n,:)'.*G(1,n).*ones(3,1);
             dNx(2,n)*eye(3) g1(n,:)'.*G(2,n).*ones(3,1)...
             g2(n,:)'.*G(2,n).*ones(3,1);
             dNx(3,n)*eye(3) g1(n,:)'.*G(3,n).*ones(3,1)...
             g2(n,:)'.*G(3,n).*ones(3,1)];
  BNL(:,ct)=A*GNL(:,ct);
end
B=BL+BNL;
epsEt=[(L(1)+ L(1)+oL(1)* L(1)+oL(4)*L(4)+oL(7)*L(7)+L(1)*oL(1)+...
        L(4)*oL(4)+ L(7)*oL(7)+ L(1)*L(1)+ L(4)*L(4)+L(7)* L(7))/2;
       (L(5)+ L(5)+oL(2)* L(2)+oL(5)*L(5)+oL(8)*L(8)+L(2)*oL(2)+...
        L(5)*oL(5)+ L(8)*oL(8)+ L(2)*L(2)+ L(5)*L(5)+L(8)* L(8))/2;
       (L(9)+ L(9)+oL(3)* L(3)+oL(6)*L(6)+oL(9)*L(9)+L(3)*oL(3)+...
        L(6)*oL(6)+ L(9)*oL(9)+ L(3)*L(3)+ L(6)*L(6)+L(9)* L(9))/2;
        L(2)+ L(4)+oL(1)* L(2)+oL(4)*L(5)+oL(7)*L(8)+L(1)*oL(2)+...
        L(4)*oL(5)+ L(7)*oL(8)+ L(1)*L(2)+ L(4)*L(5)+L(7)* L(8);
        L(8)+ L(6)+oL(3)* L(2)+oL(6)*L(5)+oL(9)*L(8)+L(3)*oL(2)+...
        L(6)*oL(5)+ L(9)*oL(8)+ L(3)*L(2)+ L(6)*L(5)+L(9)* L(8);
        L(3)+ L(7)+oL(1)* L(3)+oL(4)*L(6)+oL(7)*L(9)+L(1)*oL(3)+...
        L(4)*oL(6)+ L(7)*oL(9)+ L(1)*L(3)+ L(4)*L(6)+L(7)* L(9)];
function [N]=shapefunc(xsi,eta) 
N(:,1)= xsi.*(xsi-1).*eta.*(eta-1)     /4;
N(:,2)=-xsi.*(xsi-1).*(eta+1).*(eta-1) /2;
N(:,3)= xsi.*(xsi-1).*eta.*(eta+1)     /4;
N(:,4)=-(xsi+1).*(xsi-1).*eta.*(eta+1) /2;
N(:,5)= xsi.*(xsi+1).*eta.*(eta+1)     /4;
N(:,6)=-xsi.*(xsi+1).*(eta+1).*(eta-1) /2;
N(:,7)= xsi.*(xsi+1).*eta.*(eta-1)     /4;
N(:,8)=-(xsi+1).*(xsi-1).*eta.*(eta-1) /2;
N(:,9)=(xsi+1).*(xsi-1).*(eta+1).*(eta-1);

function [dNr]=dershapefunc(xsi,eta)
r2=size(xsi,1)*2;
dNr(1:2:r2  ,1)= eta.*(eta-1).*(2*xsi-1)    /4;
dNr(1:2:r2  ,2)=-(eta+1).*(eta-1).*(2*xsi-1)/2;
dNr(1:2:r2  ,3)= eta.*(eta+1).*(2*xsi-1)    /4;
dNr(1:2:r2  ,4)=-eta.*(eta+1).*xsi            ;
dNr(1:2:r2  ,5)= eta.*(eta+1).*(2*xsi+1)    /4;
dNr(1:2:r2  ,6)=-(eta+1).*(eta-1).*(2*xsi+1)/2;
dNr(1:2:r2  ,7)= eta.*(eta-1).*(2*xsi+1)    /4;
dNr(1:2:r2  ,8)=-eta.*(eta-1).*xsi            ;
dNr(1:2:r2  ,9)= 2*(eta+1).*(eta-1).*xsi      ;
dNr(2:2:r2+1,1)= xsi.*(xsi-1).*(2*eta-1)    /4;
dNr(2:2:r2+1,2)=-xsi.*(xsi-1).*eta            ;
dNr(2:2:r2+1,3)= xsi.*(xsi-1).*(2*eta+1)    /4;
dNr(2:2:r2+1,4)=-(xsi+1).*(xsi-1).*(2*eta+1)/2;
dNr(2:2:r2+1,5)= xsi.*(xsi+1).*(2*eta+1)    /4;
dNr(2:2:r2+1,6)=-xsi.*(xsi+1).*eta            ;
dNr(2:2:r2+1,7)= xsi.*(xsi+1).*(2*eta-1)    /4;
dNr(2:2:r2+1,8)=-(xsi+1).*(xsi-1).*(2*eta-1)/2; 
dNr(2:2:r2+1,9)= 2*(xsi+1).*(xsi-1).*eta      ;