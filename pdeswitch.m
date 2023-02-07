function pdeswitch

%---------------------------------heater model----------------------------

   h = 2;
%      h = 1;
    deltat = 0.1;
    m = 10/h + 1; %number of mesh points
    t = 10/deltat + 1; %number of time steps
    alpha = deltat/(h*h);
    A = zeros(m -2, m - 2);     %exclude dirichlet BC

    %build matrix A
    for i = 1 : m - 2
        A(i, i) = 1 - 2 * alpha;
        if i > 1
            A(i, i - 1) = alpha;
        end

        if i < m - 2
        A(i, i + 1) = alpha;
        end
    end


    disp(A)

    u0 = ones(m - 2, 1);
    u0 = u0 * 0.2;

    function f = input_fun(x)
        f = -0.1*x + 1;
    end


    load = zeros(m - 2, 1);
    u = zeros(m - 2, t);
    u(:,1) = u0;
    
    
    for i = 2 : t

        for j = 1 : m - 2	% heater on and off control

            if i == 2   %initial condition
                load(j, 1) = deltat * input_fun(h * j);
            end    

            if load(j, 1) == 0
                load(j, 1) = 0;
            else
                load(j, 1) = deltat * input_fun(h * j);
            end    
            %switching part    
            if u(j, i - 1) >= 0.7
                %u(j, i - 1) = 0.7;
                load(j, 1) = 0;
            elseif u(j, i - 1) <= 0.4
                %u(j, i - 1) = 0.4;
                load(j, 1) = deltat * input_fun(h * j);
            end
        end
        %disp(load);

        u(:, i) = A * u(:, i - 1) + load;

    end

    xlist = linspace(0,10,m);
    tlist = linspace(0,10,t);
    u = [zeros(1, t); u; zeros(1, t)];
    u = u';


%     size(xlist);
%     size(tlist);
%     size(u);



    surf(xlist,tlist,u) 
    title('Numerical solution computed with 100 mesh points.')
    xlabel('Distance x')
    ylabel('Time t')

    figure
    plot(xlist,u(91,:))
    title('Solution at t = 9')
    xlabel('Distance x')
    ylabel('u(x,5)')
    ylim([0 0.7])
    %legend('cos(x)','cos(2x)')

    figure
    plot(tlist,u(:,1 + 2/h))
    title('Solution at x = 2')
    xlabel('Time t')
    ylabel('u(0,t)')

    figure
    plot(tlist,u(:,1+4/h))
    title('Solution at x = 4')
    xlabel('Time t')
    ylabel('u(0,t)')

end