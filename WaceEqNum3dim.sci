clc; 
//функция численного решения
function [u,x,y,t] = WaveEqNum2dim(NX,NY,K,L,M,T,a)
hx=L/NX;
hy=M/NY;
tau=T/K;

for j=1:NY+1
    y(j)=(j-1)*hy;
    for i=1:NX+1
        x(i)=(i-1)*hx;
        u(i,j,1)=u_0(x(i),y(j));
        u(i,j,2)=u(i,j,1)+tau*v_0(x(i),y(j));
    end
end


for k=1:K+1
    t(k)=(k-1)*tau;
end

for j=1:NY+1
    for k=2:K+1
        u(1,j,k)=u_x0(t(k));
        u(NX+1,j,k)=u_xl(t(k));
    end
end

for i=1:NX+1
    for k=2:K+1
        u(i,1,k)=u_y0(t(k));
        u(i,NY+1,k)=u_ym(t(k));
    end
end

gx=a^2*tau^2/hx^2;
gy=a^2*tau^2/hy^2;

for k=2:K
    for j=2:NY
        for i=2:NY
            //u(i,j+1)=gam*(u(i+1,j)-2*u(i,j)+u(i-1,j))+2*u(i,j)-u(i,j-1);
            u(i,j,k+1)=gx*(u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))+gy*(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))+2*u(i,j,k)-u(i,j,k-1);
        end
    end
end
endfunction

//граничные условия
function u1 = u_x0(t)
    u1 = 0
endfunction

function u2 = u_xl(t)
    u2 = 0
endfunction

function u1 = u_y0(t)
    u1 = 0
endfunction

function u2 = u_ym(t)
    u2 = 0
endfunction

//начальные условия
function u0 = u_0(x,y)
    u0=0
endfunction

function v0 = v_0(x,y)
    if (0.4<=x)&(x<=0.6)&(0.4<=y)&(y<=0.6) then
        v0 = -2
    else
        v0 = 0
    end
endfunction


L=1;
M=1;
T=2;
a=2;
//tH=1/(a*80);
//lH=1/40;
//mH=1/40;
tH=0.01;
lH=tH*a+tH*2;
mH=tH*a+tH*2;
disp(lH,mH);
timer();
[U,X,Y,T]=WaveEqNum2dim(round(L/lH),round(M/mH),round(T/tH),L,M,T,a);
time = timer();
disp(time);
Max=max(max(U));
Min=min(min(U));

h_fig = figure;
h_fig.background = 8;
drawlater();
h_pat = plot3d(X',Y', U(:,:,1));;
h_pat.color_mode = 4;
h_pat.foreground = 1;
h_pat.hiddencolor = 4;

xlabel("X"); ylabel("Y"); zlabel("U");
h_axes = gca();
h_axes.isoview = "on";
h_axes.box = "off";
h_axes.rotation_angles = [105, -80];
h_axes.data_bounds = [0, 0, Min; 1, 1, Max];
xgrid;

outgif = 'membrana.gif';
mdelete(outgif);
idGif = animaGIF(gcf(), outgif, 10, 0);

for j=2:(length(T)*1)
    drawlater();
    p=int((j-1)./1)+1;
    h_pat.data.x = X';
    h_pat.data.y = Y';
    h_pat.data.z = U(:,:,p);
    drawnow();
    idGif = animaGIF(gcf(), idGif);
end
animaGIF(idGif);
