function r = rect(x);
  r = 0.5 * (sign(x + 0.5) - sign(x - 0.5));
end