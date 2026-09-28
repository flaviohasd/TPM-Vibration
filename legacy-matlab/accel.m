function a = accel(states,K,C,F,me)
    a(1) = states(2);
    a(2) = -K*states(1)/me -C*states(2)/me + F/me;
end