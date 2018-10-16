[coord,etpl,fext0,bc,lstps,NRitmax,NRtol,E,nu,ka]=inputfile;
nels=size(etpl,1); nodes=size(coord,1); nDoF=nodes*5;
neDoF=(9*5)^2; krow=zeros(neDoF*nels,1); kcol=krow; kval=krow;
uvw=zeros(nDoF,1); uvwold=uvw; fint=uvw; react=uvw;
fd=(1:nDoF); fd(bc(:,1))=[];
epsEn=zeros(6,8,nels); epsE=epsEn; sigN=epsEn; sig=epsEn;
Vn=zeros(9,3,nels); oVn=Vn; oL=zeros(8,9,nels); L=oL;
for lstp=0:lstps
  fext=(lstp/lstps)*fext0; foob=react+fext-fint;
  foobnorm=2*NRtol; NRit=0; 
  while ((NRit<NRitmax)&&(foobnorm>NRtol))
    NRit=NRit+1; fint=zeros(nDoF,1); dreact=fint; dduvw=fint;
    if lstp>=1
      Kt=sparse(krow,kcol,kval,nDoF,nDoF);
      dduvw(bc(:,1))=(1+sign(1-NRit))*bc(:,2)/lstps;
      dduvw(fd)=Kt(fd,fd)\(foob(fd)-Kt(fd,bc(:,1))*dduvw(bc(:,1)));
      dreact(bc(:,1))=Kt(bc(:,1),:)*dduvw-foob(bc(:,1));
    end
    uvw=uvw+dduvw; react=react+dreact; duvw=uvw-uvwold;
    for nel=1:nels 
      ed=reshape(ones(5,1)*etpl(nel,:)*5-(5-1:-1:0).'*ones(1,9),1,9*5);
      if lstp==0
        elcoord=coord(etpl(nel,:),:);
        phi=elcoord(:,4); psi=elcoord(:,5);
        oVn(:,:,nel)=[cos(psi) sin(psi).*cos(phi) sin(psi).*sin(phi)];
      end
      [ke,felem,epsE(:,:,nel),Vn(:,:,nel),sig(:,:,nel),L(:,:,nel)]=...
      shell(coord(etpl(nel,:),:),uvw(ed),duvw(ed),epsEn(:,:,nel),...
            oVn(:,:,nel),E,nu,ka,sigN(:,:,nel),oL(:,:,nel));
      if lstp==0
        ct=(nel-1)*neDoF+1:nel*neDoF;
        krow(ct)=reshape(ed.'*ones(1,9*5),neDoF,1);
        kcol(ct)=reshape(ones(9*5,1)*ed  ,neDoF,1);
      end
      kval((nel-1)*neDoF+1:nel*neDoF)=reshape(ke,neDoF,1);
      fint(ed)=fint(ed)+felem;
    end
    foob=fext+react-fint; foobnorm=norm(foob)/norm(fext+react+eps);
    fprintf('%4i %4i %6.3e\n',lstp,NRit,foobnorm);
  end
  uvwold=uvw; epsEn=epsE; sigN=sig; oL=oL+L;
end