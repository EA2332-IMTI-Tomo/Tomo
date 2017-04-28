;;(load "/users/these/bailleul/These/OpenSource/CommonLisp/Maths.sparcf")


(defun mean (values-vec
	     &key (zero 0) (add-fun #'+) (div-fun #'/))
  "return the mean of values of an (unmodified) nonempty vector/list"
  (let ((size (length values-vec))
	(sum zero))
    (assert (not (zerop size)))
    (progn (map-into values-vec 
		     #'(lambda (x) (setf sum (funcall add-fun sum x)) x)
		     values-vec)
	   (funcall div-fun sum size))))



(defun median (values-vec)
  "the x-limit cutting an histogram in two areas of equal surface"
  (let* ((sorted-vec (sort values-vec #'>))
	 (size (length sorted-vec)))
    (if (evenp size)
	(let ((vals (subseq sorted-vec (1- (/ size 2)) (1+ (/ size 2)))))
	  (/ (+ (car vals) (cadr vals)) 2))
	(let* ((middle (floor (/ size 2)))
	       (number (subseq sorted-vec middle (1+ middle))))
	  (if (typep values-vec 'cons)
	      (car number)
	      (aref number 0))))))
    
	 

(defun variance (values-vec &key (given-mean nil)
			    (zero 0) (add-fun #'+) (sub-fun #'-)
			    (sqr-fun #'(lambda (x) (* x x))) (div-fun #'/)) 
  "if mean previously computed, can be provided through given-mean to
avoid unecessary computation"
  (let ((size (length values-vec))
	(meanvalue (if (eq given-mean nil) 
		       (mean values-vec :zero zero :add-fun add-fun :div-fun div-fun)
		       given-mean))
	(cumul_sqr_dif zero))
    (assert (not (zerop size)))
    (progn (map-into values-vec
		     #'(lambda (x) 
			 (setf cumul_sqr_dif
			       (funcall add-fun cumul_sqr_dif
					(funcall sqr-fun (funcall sub-fun x meanvalue))))
			 x)
		     values-vec)
	   ;;moyenne du carre des ecarts a la moyenne
	   (funcall div-fun cumul_sqr_dif size)))) 


;;ecart-type in french
(defun standard-deviation (values-vec &key (given-mean nil))
  "return the standard deviation of given data list"
  (sqrt (variance values-vec :given-mean given-mean)))

    
(defun outliers-positions (values-vec &key (given-mean nil) (given-stddev nil))
  "return a list of indexes of outliers."
  (let* ((meanvalue (if (eq given-mean nil) (mean values-vec)
			given-mean))
	 (stddev (if (eq given-stddev nil)
		     (standard-deviation values-vec :given-mean meanvalue)
		     given-stddev))
	 (size (length values-vec))
	 (outlier-p (lambda (x) (not (boundedp x (- meanvalue stddev)
					       (+ meanvalue stddev))))) 
	 (positions-lst (list)))
    (dotimes (i size (nreverse positions-lst))
      (when (funcall outlier-p (aref values-vec i))
	(push i positions-lst)))))


(defun array-single-iterator (tab ar2fun)
  "frame to make a min/max function onto an array"
  (assert (not (zerop (length tab))))
  (let ((max (aref tab 0)))
    (map-into tab #'(lambda (x) (when (funcall ar2fun x max) (setf max x))
		      x)
	      tab)
    max))


(defun array-min (tab)
  "return minimal element of given array"
  (array-single-iterator tab #'<))


(defun array-max (tab)
  "return maximal element of given array"
  (array-single-iterator tab #'>))


(defun list-dif (lst)
  "compute the derivate of given values: return list of 1- size"
  (let ((res (list))
	(prev (pop lst)))
    (dolist (x lst (nreverse res))
      (push (- x prev) res)
      (setf prev x))))
  

(defun list-peaks (lst)
  "return a list of booleans indicating wether each list element is a
peak or not"
  (let ((len (length lst)))
    (when (endp lst) (return-from list-peaks lst))
    (when (= len 1) (return-from list-peaks '(NIL)))
    (flet ((peakp (x b a) (and (>= x b) (>= x a))))
      (let* ((last2 (lasts 2 lst))
	     (midlist (if (> len 2) (mapcar #'peakp (cdr lst) (cddr lst) lst)
			  (list)))
	     
	     (start (> (car lst) (cadr lst)))
	     (end (> (cadr last2) (car last2))))
	(cons start (pushend end midlist))))))


(defun define-gaussian (mean variance)
  "returns a gaussian function given a mean and variance"
  (flet ((square (x) (* x x)))
    (lambda (x) (exp (/ (- (square (abs (- x mean))))
			(* 2 (square variance)))))))
		      

(defun normalize (values-lst)
  "normalize given list (0/1)"
  (let* ((sum (apply #'+ values-lst)))
    (mapcar #'(lambda (x) (/ x sum)) values-lst)))


(defun percentize (values-lst)
  "normalize given list in %"
  (let* ((sum_pct (/ (apply #'+ values-lst) 100)))
    (mapcar #'(lambda (x) (/ x sum_pct)) values-lst)))


(defun percentize-cumulate (values-lst)
  "normalize given list in % and return cumulated values"
  (let* ((sum_pct (/ (apply #'+ values-lst) 100))
	 (pct-list (mapcar #'(lambda (x) (/ x sum_pct)) values-lst))
	 (cumulated (list))
	 (sum 0))
    (dolist (x pct-list (nreverse cumulated))
      (progn (incf sum x)
	     (push sum cumulated)))))


(defun amplitude (values-lst)
  (- (apply #'max values-lst)
     (apply #'min values-lst)))


;;context: put all values between 0 and 1 
(defun frame01 (values-lst)
  "bounds given values between 0 and 1"
  (let* ((minval (apply #'min values-lst))
	 (ampl (amplitude values-lst)))
    (mapcar #'(lambda (x) (no/0 ampl (/ x ampl)))
	    (mapcar #'(lambda (x) (- x	minval))
		    values-lst))))


;;context: put some focused values between 0 and 1, while discarding
;;other intermediate/extreme values that should not be accounted
;;because nonrepresentative.
;;
;;1st list: values to process, 2nd: sublist of values to account for
;;determining transformation
;;
;;important: resulting discarded values might be outside [0, 1]
(defun frame01-focused (lstvalues-toprocess lstvalues-toaccount)
  "like frame01, but values are recomputed using min and amplitude of
another list"
  (let* ((minval (apply #'min lstvalues-toaccount))
	 (ampl (amplitude lstvalues-toaccount)))
    (mapcar #'(lambda (x) (no/0 ampl (/ x ampl)))
	    (mapcar #'(lambda (x) (- x	minval))
		    lstvalues-toprocess))))


(defun list-cumul (l)
  "return a list of cumulated values of l:(a b c) -> (a a+b a+b+c)"
  (let ((res (list)))
    (dolist (e l (nreverse res))
      (push (+ e (if res (car res) 0))
	    res))))


(defun bound-values (l min max)
  "returns a list where values are bound to given min and max: if a
value is out of given bounds, it is replaced by list-min/max"
  (let* ((nolower (remove-if #'(lambda (x) (< x min)) l))
	 (noupper (remove-if #'(lambda (x) (> x max)) l))
	 (minval (if (= (length nolower) 0) min
		     (apply #'min nolower)))
	 (maxval (if (= (length noupper) 0) max
		     (apply #'max noupper))))
    (mapcar #'(lambda (x) (cond ((> x max) maxval)
				((< x min) minval)
				(t x)))
	    l)))


(defun bound01 (listvals)
  "truncates values to 0/1 if exceed [0,1]"
  (mapcar #'(lambda (x) (frame-value x 0 1))
	  listvals))


;; cf. also Maths.lisp