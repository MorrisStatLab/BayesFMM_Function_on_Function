function PI=Update_pi(gamma,a_pi,b_pi,wavespecs)
% This function update PI from beta distribution.
% Here assume PI depend on i,j, so does a_pi, b_pi.

sumop=uneqkron(wavespecs.Kj);
PI=betarnd(gamma*sumop+a_pi,(1-gamma)*sumop+b_pi);

