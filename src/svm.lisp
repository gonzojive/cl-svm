(in-package :cl-svm)

(defun dot-product (v1 v2)
  (reduce #'+ (map 'vector #'* v1 v2) :initial-value 0))

(defmacro training-pair-bind ((x-var y-var) training-pair &body body)
  "Binds the symbols named by x-var and y-var to the x and y components of a training pair."
  `(destructuring-bind (,x-var . ,y-var) ,training-pair ,@body))

(defun simple-train-svm (training-pairs &key (kkt-tolerance .005) (kernel 'dot-product) (C .1))
  "Trains a simple SVM using a simple dot-product kernel.  Training-pairs is a sequence
of of pairs, where the car of each pair is a vector of numbers and the cdr is either -1 or 1.

This will return three values:  The first value is a function that can be called with
an argument of the same dimension as a training example.  It will return either t
or nil.

The other two return values are the vector w and the scalar b that define the plane of the form
w . x - b  = 0"
  ;; The form of the optimization we wish to solve is the following dual problem:
  ;;   Max over alpha params of
  ;;     W(alpha params) = Sum(alpha) - 1/2 Sum over i,j(yi yj alphai alphaj K(xi, xj))
  ;;     s.t. 0 <= alphai <= C, Sum(alphai yi) = 0; for i = 1, 2, ..., m
  ;;
  ;; We use the SMO algorithm to perform the optimization over the alpha parameters.
  ;; The basic form of the SMO algorithm is to jointly optimize two parameters
  ;; of a constrainted optimization goal.
  ;; Repeat until convergence:
  ;;    1.  Select some pair of parameters to update next (using some heuristic so that
  ;;        we pick parameters that lead us to a global optimum)
  ;;    2.  Reoptimize W(params) with respect to the two parameters we chose while holding
  ;;        all other parameters fixed
  ;; To test for convergence while optimizing the dual problem posed by optimizing SVMs,
  ;; we can check weather the Karush-Kuhn-Tucker (KKT) conditions are satisfied within some tolerance
;  (declare (optimize debug))
  (let* ((training-pairs-count (length training-pairs))
	 ;; initialize the alpha parameters to 0.  This satisfies the 
	 (alphas (make-array training-pairs-count :initial-element 0.0))
	 (b 0.0))

    (labels ((apply-kernel (x y) (funcall kernel x y))
	     (test-convergence ()
	       (every #'(lambda (training-pair alpha)
			  (training-pair-bind (x y) training-pair
			    (let* ((w nil)
				   (val (* y (+ (apply-kernel w x) b))))
			      (cond
				((<= alpha kkt-tolerance)	  (>= val (- 1 kkt-tolerance)))
				((>= alpha (- C kkt-tolerance))   (<= val (+ 1 kkt-tolerance)))
				(t                           (and (>= val (- 1 kkt-tolerance))
								  (<= val (+ 1 kkt-tolerance))))))))
		      training-pairs alphas))
	     (random-not-n (max-val n)
	       (do ((i (random max-val) (random max-val)))
		   ((not (= i n)) i)))
	     (choose-alphas-to-optimize ()
	       "Returns two indexes into the alpha array that we should optimize"
	       (let* ((i (random training-pairs-count))
		      (j (random-not-n training-pairs-count i)))
		 (values i j)))
	     (predicted-y (x)
	       "Uses the current alphas and b to determine the predicted y value of the given vector"
	       (let ((value 0.0))
		 (map nil #'(lambda (training-pair alpha)
			      (training-pair-bind (xi yi) training-pair
				(incf value (+ (* alpha yi (apply-kernel xi x)) b))))
		      training-pairs alphas)
		 value)))

      (do ((converged nil (test-convergence)))
	  (converged
	   (values
	    #'predicted-y
	    alphas
	    b))
	(format t "Alphas: ~A~%" alphas)
	;; first we choose the alphas to jointly optimize
	(multiple-value-bind (i j) (choose-alphas-to-optimize)
	  (training-pair-bind (xi yi) (elt training-pairs i)
	    (training-pair-bind (xj yj) (elt training-pairs j)
	      (let ((alphai (elt alphas i))
		    (alphaj (elt alphas j)))
		;; next we find bounds L and H such that L <= alphaj <= H must hold
		;; in order for 0 <= alphaj <= C to hold
		(multiple-value-bind (l h)
		    (if (= yi yj)
			(values (max 0 (+ alphai alphaj (- C)))
				(min C (+ alphai alphaj)))
			(values (max 0 (- alphaj alphai))
				(min C (+ C alphaj (- alphai)))))
		  ;; compute the new value for alphaj
		  (let* ((ei (- (predicted-y xi) yi))
			 (ej (- (predicted-y xj) yj))
			 (xi.xj (apply-kernel xi xj))
			 (xi.xi (apply-kernel xi xi))
			 (xj.xj (apply-kernel xj xj))
			 (n     (+ (* 2.0 xi.xj)
				   (- xi.xi)
				   (- xj.xj)))
				 
			 (new-alphaj
			  (let* ((unclipped
				  (unless (< (abs n) .0000001)
				    (- alphaj
				       (/ (* yj (- ei ej))
					  n)))))
			    (cond
			      ((null unclipped)    (error "we cannot handle null n: ~A" n))
			      ((> unclipped h) h)
			      ((< unclipped l) l)
			      (t unclipped))))
			 (new-alphai (+ alphai (* yi yj (- alphaj new-alphaj))))
			 (new-b
			  (let ((b1
				 (- b ei (* yi (- new-alphai alphai) xi.xi) (* yj (- new-alphaj alphaj) xi.xj)))
				(b2
				 (- b ej (* yi (- new-alphai alphai) xi.xj) (* yj (- new-alphaj alphaj) xj.xj))))
			  (cond
			    ((and (< 0 new-alphai) (< new-alphai C))
			     b1)
			    ((and (< 0 new-alphaj) (< new-alphaj C))
			     b2)
			    (t
			     (/ (+ b1 b2) 2.0))))))
		    (setf b new-b
			  (elt alphas i) new-alphai
			  (elt alphas j) new-alphaj)))))))))))
			     
		    
		    
			  
		    
		    
      
    
    
	 
	 

  
  
  

