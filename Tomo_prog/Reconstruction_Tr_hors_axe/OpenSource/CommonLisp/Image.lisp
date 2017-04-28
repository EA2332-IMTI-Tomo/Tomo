;;(load "/users/these/bailleul/These/OpenSource/CommonLisp/Tools.lisp")
;;(load "/users/these/bailleul/These/OpenSource/CommonLisp/Maths.lisp")


;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;;------------------------------------------------------------------------------
;; Histogram structure
;;------------------------------------------------------------------------------

;;min not certified for now

;; an histogram is initialized from a list of colors (of course,
;; possibly non-unique)
(defclass histogram ()
  ((colors-list :initarg :colors-list :reader get-list)
   (array :accessor get-array)
   (size :accessor get-size)
   (min-color :accessor get-min)
   (max-color :accessor get-max)))


;;we allocate from 0 to max, whatever min
(defmethod initialize-instance :after ((h histogram) &key)
  (setf (get-min h) (apply #'min (get-list h)))
  (setf (get-max h) (apply #'max (get-list h)))
  (setf (get-array h)
	(let ((histogram-array (make-array (1+ (get-max h)) :element-type 'integer :initial-element 0)))
	  (progn (mapcar (lambda (i) (incf (aref histogram-array i))) (get-list h))
		 histogram-array)))
  (setf (get-size h) (length (get-array h))))


;;we allocated array from 0 to max, whatever min
(defmethod color-count ((h histogram) color)
  "retrieve color count of given color from histogram. out-of-bounds
color will raise 0 (logical, indeed)"
  (if (boundedp color (get-min h) (get-max h))
      (aref (get-array h) color) ;; (- color (get-min h))
      0))


(defmethod colors-list ((h histogram))
  "raises a list of colors inferred from histogram"
  (let ((acc (list)))
    (loop for i from (get-min h) to (get-max h)
	  do (dotimes (j (color-count h i) nil)
	       (push i acc)))
    (nreverse acc)))


;;------------------------------------------------------------------------------
;; Histogram distances
;;------------------------------------------------------------------------------


(defun distance-l2norm (histo1 histo2)
  "distance between histograms: sum of squared diffs between each color count"
  (let ((sum 0))
    (loop for i from (min (get-min histo1) (get-min histo2))
	  to (max (get-max histo1) (get-max histo2))
	  do (let ((val1 (color-count histo1 i))
		   (val2 (color-count histo2 i)))
	       (incf sum (square (- val1 val2)))))
    sum))


(defun distance-hellinger (histo1 histo2)
  "distance between histograms: sum of squared diffs between each
sqrt'ed color count"
  (let ((sum 0))
    (loop for i from (min (get-min histo1) (get-min histo2))
	  to (max (get-max histo1) (get-max histo2))
	  do (let ((val1 (color-count histo1 i))
		   (val2 (color-count histo2 i)))
	       (incf sum (square (- (sqrt val1) (sqrt val2))))))
    sum))


(defun align-histograms (histo1 histo2 distance-fun)
  "from which value shall we increase h1 colors such as (distance-fun
h1, h2) becomes minimal? For each possible value, return a cons of it
and of resulting distance-fun"
  (let* ((h1-min (get-min histo1))
	 (h2-min (get-min histo2))
	 (hd (if (<= h1-min h2-min) histo1 histo2))
	 (hf (if (<= h1-min h2-min) histo2 histo1))
	 (sign (if (eq histo1 hd) 1 -1))
	 (hd-list (colors-list hd))
	 (dist-list (list)))
    (loop for i from 0
	  to (- (get-max hf) (get-min hd))
	  do (let* ((colors-i (mapcar #'(lambda (x) (+ x i)) hd-list))
		    (hd-i (make-instance 'histogram :colors-list colors-i)))
	       ;;(break)
	       (push (cons (* i sign) (funcall distance-fun hd-i hf))
		     dist-list)))
    (nreverse dist-list)))
	  

(defun elect-histograms-alignments (conslist)
  "among possible h. a. hypothesis, raises the one with minimal
distance and (if necessary) minimal displacement"
  (car (sort conslist
	     #'(lambda (cs1 cs2)
		 (let ((val1 (cdr cs1))
		       (val2 (cdr cs2)))
		   (if (= val1 val2)
		       (< (abs (car cs1)) (abs (car cs2))) 
		       (< val1 val2)))))))


;;------------------------------------------------------------------------------
;; Binarization through mean
;;------------------------------------------------------------------------------


;; purpose: being used by binarization functions
(defun meanBW (myhistogram threshold)
  "return a pair of mean B/W values inferred from thresholding"
  (let ((black 0) (white 0) (blackhits 0) (whitehits 0))
    (loop for i from (get-min myhistogram) to (get-max myhistogram)
	  do (let ((curval (color-count myhistogram i)))
	       (when curval
		 (if (<= i threshold)
		     (progn (incf black (* curval i)) (incf blackhits curval))
		     (progn (incf white (* curval i)) (incf whitehits curval)))))
	  finally (return (cons (no/0 blackhits (/ black blackhits))
				(no/0 whitehits (/ white whitehits)))))))


;; there exists cases where t close from middle, but out of histogram bounds
(defun binarize-histogram-1 (myhistogram)
  "finds a threshold t using the following principle: for every
possible t, compute the mean B/W values, and finally select the t
closest to the middle of B and W."
  ;;no need to search for threshold outside bounds
  (let* ((tlist (iota-range (get-min myhistogram) (get-max myhistogram))) 
	 ;; for each possible t, compute B/W means 
	 (hypothesis (mapcar #'(lambda (thr)
				 (meanBW myhistogram thr))
			     tlist))	
	 ;; compute differences from mean of segment
	 (offsets (mapcar #'(lambda (bwpair thr)
			      (let* ((min (car bwpair))
				     (max (cdr bwpair))
				     (middle (+ min (/ (- max min) 2))))
				(abs (- thr middle))))
			  hypothesis tlist)))
    (+ (position (apply #'min offsets) offsets)
       (get-min myhistogram))))
  

;;------------------------------------------------------------------------------
;; Binarization through mean - Jol
;;------------------------------------------------------------------------------


(defun threshold-start (array-histo &optional (min 0) (max 255))
  (declare (type '(simple-array fixnum 256) array-histo))
  (let ((left min) (right max) (total 0) (i 0))
    (declare (type fixnum left right total i))
    (while (zerop (aref array-histo i))
      (incf i))
    (setf left i
          i max)
    (while (zerop (aref array-histo i))
      (decf i))
    (setf right i)
    (loop for ii from left to right do
	  (incf total (aref array-histo ii)))
    (list total left right)))


(defun binarize-histogram-2 (histogram-obj &optional (min 0) (max (1- (get-size histogram-obj))))
  (let* ((myhistogram (get-array histogram-obj))
	 (ts (threshold-start myhistogram min max))
	 (total (first ts))
         (left (second ts))
         (right (third ts))
         (hl (aref myhistogram left))
         (bl hl)			; bl(letf)
         (br (- total hl))		; br(left)
         (al (/ hl 2))			; al(left)
         (ar 0)
         (best-threshold (1+ left))
         (sum-from-0 hl)
         best-value
         hi)
    ;;(scores (make-array 256 :element-type 'float :initial-element 0.0)))
    (loop for i from (1+ left) to right do
          (incf ar (* (- i left 0.5) (aref myhistogram i)))) ; ar(left)
    ;;(setf best-value (abs (- (/ al bl) (/ ar br))))
    ;;(setf (aref scores (1+ left)) best-value)
    (setf best-value (* (/ al bl) (/ ar br)))
    ;;(with-open-file (vals "vals.txt" :direction :output :if-exists :supersede)
    ;;(format vals "~A ~A~%" (1+ left) best-value)
    (loop for i from (1+ left) to (1- right) do
          (setf hi (aref myhistogram i)
                al (+ al 
                      (* 0.5 hi)
                      sum-from-0))
          (incf sum-from-0 hi)
          (setf ar (- ar (* 0.5 hi) (- total sum-from-0)))
          (incf bl hi)
          (decf br hi)
          ;;(let ((new-value (abs (- (/ al bl) (/ ar br)))))
          (let ((new-value (* (/ al bl) (/ ar br))))
            ;;(setf (aref scores (1+ i)) new-value)
            ;;(format vals "~A ~A~%" (1+ i) new-value)
            (when (> new-value best-value)
              (setf best-value new-value
                    best-threshold (1+ i))))
          finally (return (values best-threshold br))))) 



;;------------------------------------------------------------------------------
;; for easier use
;;------------------------------------------------------------------------------


(defun threshold-listdata (l &key (method 1))
  "returns the threshold of any list data"
  (let* ((data256 (mapcar #'(lambda (x) (round (* 256 x)))
			  (frame01 l)))
	 (histo (make-instance 'histogram :colors-list data256)))
    (/ (case method
	 (1 (binarize-histogram-1 histo))
	 (2 (binarize-histogram-1 histo))
	 (otherwise (error "undefined binarization method")))
       256)))
	


;;------------------------------------------------------------------------------
;; deprecated 
;;------------------------------------------------------------------------------


(defun make-array-histogram (values-lst &key (maxval nil)) 
  "makes a simplistic histogram (array) from a list of integer
values. If no max-value specified, takes the maximum of the list"
  (let* ((maxint (if maxval maxval (apply #'max values-lst)))
	 (myhistogram (make-array (1+ maxint) :element-type 'integer :initial-element 0)))
    (progn (mapcar (lambda (i) (incf (aref myhistogram i))) values-lst)
	   myhistogram)))


;;------------------------------------------------------------------------------
;; RGB
;;------------------------------------------------------------------------------

(deftype rgb () 'vec3)


(defun rgb-make (r g b)
  (vec3-normalize (vec3-makeg r g b)))
(defmacro rgb-r (c) `(vec-x ,c))
(defmacro rgb-g (c) `(vec-y ,c))
(defmacro rgb-b (c) `(vec-z ,c))


;;HSV funs
(defun rgb-saturation (c)
  (let ((cmax (vec3-max c))
	(cmin (vec3-min c)))
    (/ (- cmax cmin) cmax)))

(defmacro rgb-value (c) `(vec3-max ,c))

(defun rgb-hue (c)
  "colors in trigo order: Red, Yellow, Green, Cyan, Blue, Magenta"
  (let* ((cmax (vec3-max c))
	 (cmin (vec3-min c))
	 (r (rgb-r c))
	 (g (rgb-g c))
	 (b (rgb-b c))
	 (delta (- cmax cmin))
	 (dr (/ r delta))
	 (dg (/ g delta))
	 (db (/ b delta))
	 (temphue (* 60 (cond ((= cmax r) (- db dg)) ; r dominant
			      ((= cmax g) (+ 2 (- dr db)))
			      ((= cmax b) (+ 4 (- dg dr)))
			      (t (error "§§§"))))))
    (if (< temphue 0)
	(+ temphue 360)
	temphue)))