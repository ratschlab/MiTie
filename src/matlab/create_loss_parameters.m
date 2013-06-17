function fn_save = create_loss_parameters(eta1, eta2, lambda, dir_name)


% mean and variance for the poison noise
% lambda
%
% var = eta1*mu+eta2*mu^2

fn_save = sprintf('%s/loss_param_lambda-%i_eta1-%.2f_eta2-%.2f.mat', dir_name, lambda, eta1, eta2);

if fexist(fn_save)
	return
end

do_plot = 1;

if do_plot
	figure
end

% this is the observed coverage:
%xpos=[1 2 3 5 10 15 20 35 50 100 200 500 1000 5000 10000 30000] ;
xpos = round(logspace(0, 4, 20));

s=0; 
for obs=xpos ;
	s=s+1 ;

	std_=sqrt(eta1*(obs+1) + eta2*obs^2);
	%mus=max(1, obs-5*std_):step:(obs+5*std_);
	mus = [];
	for j = -5:4
		if  obs+j*std_<1 && j<-1
			continue
		end
		mus = [mus linspace(max(1, obs+j*std_), max(1, obs+(j+1)*std_), 20)];
	end
	
	Y4 = zeros(1, length(mus));
	Y = zeros(1, length(mus));
	for k=1:length(mus),
		mu = mus(k);
		var_ = eta1*mu + eta2*mu^2;
		[r p] = compute_rp_neg_binom(mu, var_); 
		% A ~ Pois(lambda)
		% B ~ NB(r,p)
		% obs = A+B 
		%
		for i = 0:min(obs, 100)
			if lambda==0 && i==0
				pois = 1;
			elseif lambda==0
				break;
			else
				pois = i*log(lambda) - lambda - factln(i); 
			end
			if eta1==1 && eta2==0
				% special case: in the limit of var->mu we get the poisson distrib
				nb = (obs-i)*log(mu) - mu - factln(obs-i);
			else
				nb = factln(obs-i+r-1) - factln(obs-i) - factln(r-1) + r*log(1-p) + (obs-i)*log(p); 
			end
			%Y4(k) = Y4(k) + exp(pois+nb); 

			% using log(a+b) = log(a) + log(1+exp(log(b)-log(a)))
			if i==0
				Y(k) = nb+pois;
			else
				Y(k) = Y(k) + log(1+exp(nb+pois-Y(k)));
			end
		end
	end ;
	%Y4(Y4<1e-200) = 1e-200;
	%Y = -log(Y4);
	Y = -Y;

	%clf
	%plot(mus, Y), hold on, plot(mus, -log(Y4), 'r'); 
	%keyboard
	%continue


	X=[(mus-obs).^2; abs(mus-obs)] ;
	idx1=find(mus<=obs) ;
	idx2=find(mus>=obs) ;
	%Y=-log(P') ;
	Y = Y';
	if ~isempty(idx1)
		offset=Y(idx1(end)) ;
	else
		offset=Y(idx2(1));
	end

	% fit the function locally: fit is more exact close to the 
	% observed value than far away from it
	%WWA=exp(-(mus-obs).^2./median(((mus-obs)).^2)) ;
	opts =  optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
	
	if 1,
	    %global XX YY WW
	    XX=X(:,idx1)' ;
	    YY=Y(idx1)-offset;
	    %WW=WWA(idx1) ;
		WW = ones(1, length(idx1));
		if 0%(exist('fmincon')==2)
			w1 = [0;0];
	   		w1 = fmincon(@(w) mean(WW'.*((XX*w-YY).^2)),w1,[],[], [], [], [1e-10 10], [10 1e-10]) 
		elseif isempty(idx1) || all(YY==0)
			w1 = [1e-10; 1e-10];
		else
			w1 = my_min(XX, YY, WW, [1e-10 10], [10 1e-10])
		end
	    E_left=mean(WW'.*((XX*w1-YY).^2))
	
		if 0% all(w1==0)
			figure, hold on
			plot(mus(idx1), WW), plot(mus(idx1), YY), plot(mus(idx1), XX*[1e-7;1e-1], 'g')
			keyboard
		end
	    
	    %global XX YY
	    XX=X(:,idx2)' ;
	    YY=Y(idx2)-offset ;
	    %WW=WWA(idx2);
		WW = ones(1, length(idx2));
		if 0%(exist('fmincon')==2)
			w2 = [0;0];
	    	w2 = fmincon(@(w) mean(WW'.*((XX*w-YY).^2)),w2,[],[], [], [], [1e-10 0.5], [10 1e-10])
		else
			w2 = my_min(XX, YY, WW, [1e-10 0.5], [10 1e-10])
		end
	    E_right=mean(WW'.*((XX*w2-YY).^2))
	end ;
	
	left_q(s)=w1(1);
	left_l(s)=w1(2);
	right_q(s)=w2(1);
	right_l(s)=w2(2);

	if isnan(w1(1)) || isnan(w2(1))
		keyboard
	end
	
	if do_plot
		subplot(4,5,s) ;
		%figure
		semilogy(mus, Y, mus(idx1), (X(:,idx1)'*w1)'+offset, mus(idx2), (X(:,idx2)'*w2)'+offset)
		%title(sprintf('obs=%1.2f, E_{left}=%1.4f, E_{right}=%1.4f', obs, E_left, E_right)) ;
		%plot(mus,-log(P), mus, offset+(log(mus)-log(obs)).^2) ;
	end

end ;

fprintf('observation\tleft_l\tleft_q\tright_l\tright_q\n');
for j = 1:length(xpos), 
	fprintf('%i\t%.11f\t%.11f\t%.11f\t%.11f\n', xpos(j), left_l(j), left_q(j), right_l(j), right_q(j));
end

save(fn_save, 'xpos', 'left_l', 'left_q', 'right_l', 'right_q')

return

function w_best = my_min(XX, YY, WW, box1, box2)

	eps = 1e-15;
	w = [box1(1); box2(1)];
	%best = mean(WW'.*((XX*w-YY).^2));
	n = 1/sum(exp(-XX*w)); % normalize
	best = -sum(sqrt(exp(-YY)*1/sum(exp(-YY)).*(exp(-XX*w)*n))); % negative Bhattacharyya coefficient
	best_orig = best;
	w_best = w;
	iter = 0;
	% weighted least squares
	while abs(box1(2)-box1(1))>eps || abs(box2(2)-box2(1))>eps 
		steps1 = box1(1):(box1(2)-box1(1))/10:box1(2);
		steps2 =  box2(1):(box2(2)-box2(1))/10:box2(2);

		best = best_orig;
		for val1 = steps1
			for val2 = steps2
				w = [val1; val2];
				%val = mean(WW'.*((XX*w-YY).^2));
				n = 1/sum(exp(-XX*w)); % normalize
				val = -sum(sqrt(exp(-YY)*1/sum(exp(-YY)).*(exp(-XX*w)*n)));
				if val<best
					best = val;
					w_best = w;
				end
			end
		end	
		if abs(box1(2)-box1(1))>eps
			idx = find(w_best(1)==steps1, 1, 'first');
			if idx>1 && idx<length(steps1)
				box1 = [steps1(idx-1), steps1(idx+1)];
			elseif idx>1
				box1 = [steps1(end-1), steps1(end)];
			else
				box1 = [steps1(1), steps1(2)];
			end
		end
		box2_tmp = box2;
		if abs(box2(2)-box2(1))>eps
			idx = find(w_best(2)==steps2, 1, 'first');
			if idx>1 && idx<length(steps2)
				box2 = [steps2(idx-1), steps2(idx+1)];
			elseif idx>1
				box2 = [steps2(end-1), steps2(end)];
			else
				box2 = [steps2(1), steps2(2)];
			end
		end
		iter = iter+1;
		if iter>200 % may happen for numerical reasons if eps too small
			break
		end
	end
	%exp(-YY)'*1/sum(exp(-YY)), exp(-XX*w)'*1/sum(exp(-XX*w))
	iter
return


