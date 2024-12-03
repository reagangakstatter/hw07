% Author: Reagan Gakstatter / reg0052@auburn.edu
% Date: 2024-12-11
% Assignment Name: hw07

classdef hw07
    
    methods (Static)
        function y = p1(func, y0, tspan, n_steps, method)
            % Solves the ODE y' = f(t, y) with initial condition y(t0) = y0 using the specified method (euler, rk4, midpoint).
            % over the interval tspan=[a, b] with n_steps. The function f(t, y) is provided as a function handle. 
            % 
            %:param func: function handle f(t, y) that defines the ODE y' = f(t, y)
            %:param y0: initial condition y(t0) = y0
            %:param tspan: interval [a, b] over which to solve the ODE
            %:param n_steps: number of steps to take to solve the ODE, interval size = (b-a)/n_steps.
            %:param method: string that specifies the method to use. It can be 'euler', 'midpoint', or 'rk4'
            %
            %:return: none, but plots the solution y(t) over the interval tspan
            
            % Your implementation here. Euler method is implemented for an example. Implement the other methods.

            t0 = tspan(1); 
            tf = tspan(2); 
            
            h = (tf - t0) / n_steps; 
            t = t0:h:tf; 
            y = zeros(1, length(t)); 
            
            y(1) = y0;

            for i = 1:n_steps
                switch method
                    case 'euler'
                        k1 = func(t(i), y(i));
                        y(i+1) = y(i) + h * k1;
                    case 'midpoint'
                        k1 = func(t(i), y(i));
                        k2 = func(t(i) + h/2, y(i) + h/2 * k1);
                        y(i+1) = y(i) + h * k2;
                    case 'rk4'
                        k1 = func(t(i), y(i));
                        k2 = func(t(i) + h/2, y(i) + h/2 * k1);
                        k3 = func(t(i) + h/2, y(i) + h/2 * k2);
                        k4 = func(t(i) + h, y(i) + h * k3);
                        y(i+1) = y(i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
                    otherwise
                        error('Invalid method. Choose "euler", "midpoint", or "rk4".');
                end
            end
        end

        function p2(method)
            % Test the implemented methods on the ODE
            % y' = t(y - t sin(t)) with initial condition y(0) = 1 over the interval [0, 1], [0, 3], and [0, 5] with variable step lengths. 

            % 
            % Plot the solution y(t) over the interval for various step sizes h. And plot the 
            % exact solution y(t) = t sin(t) + cos(t) over the same interval.
            %
            % Run the commands below to test the implemented methods:
            %
            %> hw07.p2('euler');
            %> hw07.p2('midpoint');
            %> hw07.p2('rk4');
            %
            % Observe the solution and error plots for the numerical solutions with different step sizes. Write your observations in the comments. 

            % Your comment here (e.g, how does the error change with step size and the time span, etc.): 
            % 
            % For Euler's method, the error decreased linearly, therefore h.
            % Error accumulated the fastest due to its low-order and grew significantly big for
            % large h values.
            % Midpoint method decreased quadratically with h. Grew slower than Euler's but would
            % still accumulate at large h's.
            % RK4 error decreased steeply, h^4. It was the most efficient for achieving better
            % accuracy at larger h values compared to Euler's and midpoint. Its error grew the
            % slowest because of its higher-order accuracy.
            % 
            
             f = @(t, y) t * (y - t * sin(t));

            tf_values = [1,3,5];

            figure('Position', [0 0 1200 2000]);
            for tf_index = 1:length(tf_values)

                t0 = 0; tf = tf_values(tf_index); y0 = 1;
                exact_sol = @(t) t .* sin(t) + cos(t);
                error = zeros(1, 8);
                h = [1e-1, 1e-1 * 10/11, 1e-1*4/5, 1e-1*0.72, 1e-1*0.64, 1e-1*1/2, 1e-1*2/5, 1e-1*0.32];
                c = {'g--', 'g-.', 'b--', 'b-.', 'm--', 'm-.', 'k--', 'k-.'};

                subplot(length(tf_values),2, tf_index * 2 - 1);
                title(['Numerical', ' solution using ', method, ' method and Exact Solutions', ' on interval [', num2str(t0), ', ', num2str(tf), ']']);
                for i = 1:length(h)
                    n_steps = (tf - t0) / h(i);
                    y=hw07.p1(f, y0, [t0, tf], n_steps, method);
                    error(i) = max(abs(y - exact_sol(t0:h(i):tf)));
                    hold on;
                    plot(t0:h(i):tf, y, sprintf('%s', c{i}), 'DisplayName', ['h = ', num2str(h(i))]);
                end

                plot(t0:h(end):tf, exact_sol(t0:h(end):tf), 'r-', 'DisplayName', 'Exact Solution');
                hold off; legend("Location", 'best'); grid on; xlabel('t'); ylabel('y(t)')

                subplot(length(tf_values),2, tf_index * 2);
                    loglog(h, error, 'b-o', 'DisplayName', 'Max Error vs. Step Size');
                    hold on;
                    loglog(h, h.^4 * error(1)/h(1)^4, 'r--', 'DisplayName', '4th order convergence');
                    loglog(h, h.^3 * error(1)/h(1)^3, 'g--', 'DisplayName', '3rd order convergence');
                    loglog(h, h.^2 * error(1)/h(1)^2, 'm--', 'DisplayName', '2nd order convergence');
                    loglog(h, h.^1 * error(1)/h(1)^1, 'k--', 'DisplayName', '1st order convergence');
                    hold off;
                    title(['Error vs. Step Size for y'' = t(y - t sin(t))', ' on [', num2str(t0), ', ', num2str(tf), ']']); xlabel('Step Size (h)'); ylabel('Error'); grid on; legend("Location", 'best');
           end
       end
       
     
    end
end