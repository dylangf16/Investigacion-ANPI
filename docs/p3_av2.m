function p3_av2()
  clc; clear; close all;
  warning('off', 'all');

  pkg load symbolic

  tol = 1e-5;
  iterMax = 1000;

  % --- Primer Sistema ---
  x_0 = [3, 2];
  f = {'x1^2 - 2*x1 - x2 + 0.5', 'x1^2 + 4*x2^2 - 4'};
  x = {'x1', 'x2'};
  [x_k, i, er_k] = newton_raphson(x_0, f, x, tol, iterMax);

  fprintf('--- Primer Sistema ---\n');
  fprintf('Se obtuvieron las aproximaciones: x1 = %f y x2 = %f\n', x_k(1), x_k(2));
  fprintf('Convergio en %d iteraciones\n', i);
  fprintf('Con un error final de: %.20f\n', er_k(end));
  fprintf('--------------------------------------------------------------------------\n');

  figure;
  loglog(1:length(er_k), er_k);
  xlabel('Iteración');
  ylabel('Error (norma de F)');
  title('Convergencia de Newton-Raphson para Primer Sistema');
  grid on;

  % --- Segundo Sistema ---
  x_0 = [1.2, -1.5];
  f = {'sin(x1) + x2*cos(x1)', 'x1 - x2'};
  x = {'x1', 'x2'};
  [x_k, i, er_k] = newton_raphson(x_0, f, x, tol, iterMax);

  fprintf('--- Segundo Sistema ---\n');
  fprintf('Se obtuvieron las aproximaciones: x1 = %f y x2 = %f\n', x_k(1), x_k(2));
  fprintf('Convergio en %d iteraciones\n', i);
  fprintf('Con un error final de: %.20f\n', er_k(end));
  fprintf('--------------------------------------------------------------------------\n');

  figure;
  loglog(1:length(er_k), er_k);
  xlabel('Iteración');
  ylabel('Error (norma de F)');
  title('Convergencia de Newton-Raphson para Segundo Sistema');
  grid on;

  % --- Tercer Sistema ---
  x_0 = [-1, -1, -1, -1];
  f = {'x2*x3 + x4*(x2+x3)', 'x1*x3 + x4*(x1+x3)', ...
       'x1*x2 + x4*(x1+x2)', 'x1*x2 + x1*x3 + x2*x3 - 1'};
  x = {'x1', 'x2', 'x3', 'x4'};
  [x_k, i, er_k] = newton_raphson(x_0, f, x, tol, iterMax);

  fprintf('--- Tercer Sistema ---\n');
  fprintf('Se obtuvieron las aproximaciones: x1 = %f y x2 = %f\n', x_k(1), x_k(2));
  fprintf('Se obtuvieron las aproximaciones: x3 = %f y x4 = %f\n', x_k(3), x_k(4));
  fprintf('Convergio en %d iteraciones\n', i);
  fprintf('Con un error final de: %.20f\n', er_k(end));

  figure;
  loglog(1:length(er_k), er_k);
  xlabel('Iteración');
  ylabel('Error (norma de F)');
  title('Convergencia de Newton-Raphson para Tercer Sistema');
  grid on;

end

function [x_k, i, er_k] = newton_raphson(x_0, f, x, tol, iterMax)
  n = length(f);

  % --- Paso 1: convertir las variables a simbólico ---
  x_sym = cell(1, n);
  for j = 1:n
    x_sym{j} = sym(x{j});
  end

  % --- Paso 2: calcular el Jacobiano simbólico ---
  J_sym = sym(zeros(n, n));
  for k = 1:n
    expr = sym(f{k});
    for j = 1:n
      J_sym(k, j) = diff(expr, x_sym{j});
    end
  end

  x_k = x_0(:);  % Asegurar vector columna
  er_k = [];

  % --- Paso 3: aproximar valores ---
  for i = 1:iterMax

    % Evaluar el vector F en x_k
    f_eval = zeros(n, 1);
    for k = 1:n
      expr = sym(f{k});
      for j = 1:n
        expr = subs(expr, x_sym{j}, x_k(j));
      end
      f_eval(k) = double(expr);
    end

    % Registrar el error actual
    er_k(end+1) = norm(f_eval, 2);

    % Criterio de parada por tolerancia
    if er_k(end) < tol
      break;
    end

    % Evaluar el Jacobiano numéricamente en x_k
    J_eval = zeros(n, n);
    for fila = 1:n
      for columna = 1:n
        entry = J_sym(fila, columna);
        for j = 1:n
          entry = subs(entry, x_sym{j}, x_k(j));
        end
        J_eval(fila, columna) = double(entry);
      end
    end

    % Resolver el sistema lineal
    y = mldivide(J_eval, f_eval);

    % Actualizar la aproximación
    x_k = x_k - y;

  end
end
