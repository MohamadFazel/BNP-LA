function[t_p, sigma_p] = IRF_fit(signals , resolution , T_min , T_max )

t=T_min:resolution:T_max;

for n=1:length(signals)
    
    IRF=signals{n}*resolution;
    IRF=IRF(IRF<T_max);
    IRF=IRF(IRF>=T_min);

    
    [Gex,~]=hist(IRF,t);


    frew=@(x) fcspddeb( [x(1), x(2) , x(3) ] , t , Gex);
    
    options = optimoptions('fmincon','TolX',1e-8,'MaxFunEvals',2000,'MaxIter',2000,'StepTolerance',1e-6);
    
    x0      = rand(1,3)                      ;
    A       = []                             ; 
    b       = []                             ; 
    Aeq     = []                             ; 
    beq     = []                             ;
    lb      = [0,T_min,0]                      ;
    ub      = [10^7,T_max,10]              ;
    nonlcon = []                             ;
    
    [x]=fmincon(frew,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    
    fcspddeb( [x(1),x(2),x(3)] , t,Gex);

    T_p(n)   = x(2);
    Sig_p(n) = x(3);

end

t_p     = mean(T_p)   ;
sigma_p = mean(Sig_p) ;

end


function [ X2 ] = fcspddeb( x,  t , Gex  )

             G  = x(1).*exp(-((t-x(2)).^2)./(2*(x(3).^2))) ;
           
             X2 = sum((Gex-G).^2);
end
