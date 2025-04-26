clc; 
//функция численного решения волнового уравнения в 2-х измерениях
function [u,x,t] = WaveEqNum(N,K,L,T,a)
h=L/N;
tau=T/K;


for i=1:N+1
    x(i)=(i-1)*h;
    u(i,1)=u_0(x(i));
    u(i,2)=u(i,1)+tau*v_0(x(i));
end


for j=1:K+1
    t(j)=(j-1)*tau;
end
for j=2:K+1
    u(1,j)=u_1(t(j));
    u(N+1,j)=u_2(t(j));
end

//gam=a^2*tau^2/h^2;

for j=2:K
    for i=2:N
        //u(i,j+1)=gam*(u(i+1,j)-2*u(i,j)+u(i-1,j))+2*u(i,j)-u(i,j-1);
        u(i,j+1)=gam(x(i),a,tau,h)*(u(i+1,j)-2*u(i,j)+u(i-1,j))+2*u(i,j)-u(i,j-1);
    end
end
endfunction

function r = c(x)
    if (0<=x) & (x<=0.5) then
        r = 1
    else
        r = 1
    end
endfunction

function g =gam(x,a,t,h)
    g = c(x)^2*a^2*t^2/h^2
endfunction
//граничные условия
function u1 = u_1(t)
    u1 = 0
endfunction

function u2 = u_2(t)
    u2 = 0
endfunction

//начальные условия
function u0 = u_0(x) //начальная форма волны в каждой точке
    //u0 = x*(x^2-1);
    //u0 = sin(x);
    /*
    if (0<=x) & (x<=5) then
        u0 = 0.6*x
    else
        u0 = -0.6*x+6
    end
    */
    //u0 = x*(x^2-1)
    u0=0
endfunction

function v0 = v_0(x) //начальная скорость волны в каждой точке
    /*
    if (0<=x) & (x<=0.5) then
       v0 = x^2-0.5*x
    else
        v0 = -x^2+1.5*x-0.5
    end
    */
    if (0<=x) & (x<=0.5) then
        v0 = 2*sin(2*%pi*x)
    else
        v0 = -0.5*sin(2*%pi*x)
    end
endfunction


L=1;
T=60;
a=0.05;
lH=0.02;
tH=lH/a;
[U,X,T]=WaveEqNum(round(L/lH),round(T/tH),L,T,a);

surf(X,T,U');
xlabel('X');
ylabel('T');
zlabel('U');


Max=max(max(U));
Min=min(min(U));
h_fig = figure; 
h_fig.background = 8;
h_point = plot(X, U(:,1));
xgrid();
gca().data_bounds = [0, L, Min, Max];

outgif = 'wave.gif';
mdelete(outgif);
idGif = animaGIF(gcf(), outgif, 10, 0);

for j=2:(length(T+1)*1)
    drawlater(); 
    p=int((j-1)./1)+1;
    h_point.data = [X,U(:,p)];
    drawnow();
    idGif = animaGIF(gcf(), idGif);
end
animaGIF(idGif);
clear;

