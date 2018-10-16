function [coord,etpl,fext,bc,lstps,NRitmax,NRtol,E,nu,ka]=inputfile
NRitmax=50; NRtol=1e-6; ndim=5; lstps=20;
nels=1; len=10; dep=0.1; wid=1;
E=1.2e6; nu=0; ka=5/6;
I=(wid*dep^3)/12; M=2*pi*E*I/len;

coord=zeros((nels*2+1)*3,6);
etpl=zeros(nels,9);
coord(:,4:5)=pi/2; 
coord(:,6)=dep;
for i=1:(nels*2+1)
  coord((i-1)*3+1,1:2)=[(i-1)*len/(2*nels) 0];  
  coord((i-1)*3+2,1:2)=[(i-1)*len/(2*nels) wid/2];
  coord((i-1)*3+3,1:2)=[(i-1)*len/(2*nels) wid];
end
for i=1:nels
  n1=(i-1)*6+1;
  n2=n1+1;
  n3=n1+2;
  n4=n1+5;
  n5=n1+8;
  n6=n1+7;
  n7=n1+6; 
  n8=n1+3;
  n9=n1+4;
  etpl(i,:)=[n1 n2 n3 n4 n5 n6 n7 n8 n9];    
end
nodes=size(coord,1);
fext=zeros(nodes*ndim,1);
fext( nodes*5   -2)=4*1/6;
fext((nodes-1)*5-2)=4*2/3;
fext((nodes-2)*5-2)=4*1/6;
% fext( nodes*5   -4)=4e4*1/6;
% fext((nodes-1)*5-4)=4e4*2/3;
% fext((nodes-2)*5-4)=4e4*1/6;
% fext( nodes*5   -0)=-50*pi/3*1/6;
% fext((nodes-1)*5-0)=-50*pi/3*2/3;
% fext((nodes-2)*5-0)=-50*pi/3*1/6;
% fext( nodes*5   -0)=-M*1/6;
% fext((nodes-1)*5-0)=-M*2/3;
% fext((nodes-2)*5-0)=-M*1/6;
bc=zeros(nodes*5,2); k=0;
for i=1:nodes
  if coord(i,1)==0
    bc(i*5-4,:)=[i*5-4 0]; k=k+1;
    bc(i*5-3,:)=[i*5-3 0]; k=k+1; 
    bc(i*5-2,:)=[i*5-2 0]; k=k+1;
    bc(i*5-1,:)=[i*5-1 0]; k=k+1;
    bc(i*5  ,:)=[i*5   0]; k=k+1;
  end
end

bc=sortrows(bc,1);
bc(1:nodes*5-k,:)=[];