function ph = p4Phases(n)
%P4PHASES P4 code (Lewis and Kretschmer): phase_k = pi*k^2/N - pi*k, k = 0..N-1.
k = (0:n-1).';
ph = pi * k.^2 / n - pi * k;
end
