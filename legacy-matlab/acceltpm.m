function a = acceltpm(states,K,C,F,me,d1,d2,K1,K2)
    a(1) = states(2);
    a(2) = (1/me)*(-K*states(1) -C*states(2) + F + K1*d1 + K2*d2);
end