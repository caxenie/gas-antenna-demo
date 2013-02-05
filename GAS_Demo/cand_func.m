function f = cand_func(x)
  f(:,1) = x(:).*sin(10*pi*x(:)) + 1;
end
