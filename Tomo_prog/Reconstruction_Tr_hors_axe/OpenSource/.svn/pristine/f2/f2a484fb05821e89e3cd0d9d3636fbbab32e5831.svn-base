;;a stack structure of unlimited size providing direct access to elements
(defclass queue ()
  ((cluster-size :accessor get-cluster-size
		 :initarg :cluster-size
		 :initform 10)
   (vector :accessor private-vector)))


(defmethod initialize-instance :after ((q queue) &key (verbose nil))
  :documentation "called after queue object instanciation"
  (progn (setf (private-vector q)
	       (make-array (get-cluster-size q) :fill-pointer 0 :initial-element nil 
			   :adjustable t))
	 (when verbose
	   (format t "~%queue created with cluster size ~A" (get-cluster-size q)))))


(defmethod q-emptyp ((q queue))
  :documentation ""
  (not (> (length (private-vector q))
	  0)))


(defmethod q-length ((q queue))
  :documentation ""
  (length (private-vector q)))


(defmethod q-push (element (q queue))
  :documentation "returns new length of queue"
  (let ((size (q-length q))
	(cluster (get-cluster-size q)))
    ;; notice test condition always evaluated
    (when (eql nil (vector-push element (private-vector q))) 
      (progn (adjust-array (private-vector q) (+ size cluster) :initial-element nil)
	     (assert (vector-push element (private-vector q)))))
    (q-length q)))


(defmethod q-push-seq (sequence (q queue))
  :documentation "returns new length of queue"
  (progn
    (map 'list
	 #'(lambda (elem) (q-push elem q))
	 sequence)
     (q-length q)))


(defmethod q-pop (element (q queue))
  :documentation "dequeues and return last element"
  (vector-pop (private-vector q)))


(defmethod q-element ((q queue) index)
  :documentation "return i-th element of queue"
  (let ((size (q-length q)))
    (assert (< index size))
    (aref (private-vector q) index)))


(defmethod q-mapcar (arity-1-fun (q queue) &key (queue-output nil))
  "applies given function to queue elements without modification using
mapcar scheme. (return a list)"
  (let ((result-list (map 'list
			  arity-1-fun
			  (private-vector q)))
	(result-queue (make-instance 'queue)))
    (if queue-output
	(progn (q-push-seq result-list result-queue)
	       result-queue)
	result-list)))
	  

(defmethod q-map-into (arity-1-fun (q queue))
  :documentation "maps given function to queue elements with queue modification (x -> f(x))."
  (map-into (private-vector q)
	    arity-1-fun
	    (private-vector q)))


(defmethod q-display ((q queue) &key (format_str "~A "))
  :documentation "displays every elements of given queue"
  (progn (q-mapcar #'(lambda (elt) (format t format_str elt))
		   q)
	 t))


(defmethod q-save ((q queue) path &key (stringize-fun #'mkstr))
  :documentation "saves contents of given queue to file. One element
per line, whose string representation computed by stringize-fun."
  (with-open-file (stream path :direction :output :if-exists :supersede)
    (map-into (private-vector q)
	      #'(lambda (vec)
		  (progn (write-line (funcall stringize-fun vec) stream)
			 vec))
	      (private-vector q)))
  t)


(defmethod q-list ((q queue))
  :documentation "return equivalent list"
  (map 'list #'identity (private-vector q)))


;;(setf q (make-instance 'queue))

;;(dotimes (i 25)
;;  (q-push i q))
;;NIL

;;(q-length q)
;;25

 