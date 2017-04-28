;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)



;;reads a vtk file represention of one 3D object and returns a
;;vec3-vector of 3d vertices contained inside.
;;;;files must be in unix form: no C-m or C-z
(defun vtk-parse-vertices (path
			   &key (only-nb_vertex nil)) ;; set to t if you only want to retrieve number of vertices
  "return vertices vector from vtk file"
  (let ((get-headtoken-secure (lambda (line) (if (equalp 'eof line) 'eof
						 (str-first-token line))))
	(vec3-vector nil)
	(vec3-vector-mark 0)
	(parsing-complete nil))
    (with-open-file
	(stream path :direction :input)
      (do* ((current-line (read-line stream nil 'eof)
			  (read-line stream nil 'eof)) ;current line of the stream
	    (first-token (funcall get-headtoken-secure current-line)
			 (funcall get-headtoken-secure current-line))
	    (state 0)
	    (nbpoints nil))
	  ((or (equal current-line 'eof) parsing-complete) vec3-vector) 
	(cond ((equalp first-token "#") 'comment)
	      ;; read introductory lines	      
	      ((and (eq state 0) (equalp first-token "vtk"))
	       (incf state))
	      ((and (eq state 1) (equalp first-token "ASCII"))
	       (incf state))
	      ((and (eq state 2) (equalp first-token "DATASET"))
	       (if (not (equalp (str-fetch-arg "DATASET" current-line) "POLYDATA"))
		   (error "vtk file parsing failed: POLYDATA only")
		   (incf state)))
	      ;; nb points read: allocate vec3-vector	      
	      ((and (eq state 3) (equalp first-token "POINTS"))
	       (setf nbpoints (str-coerce-i (str-fetch-arg "POINTS" current-line)))
	       (when only-nb_vertex (return nbpoints))
	       (format t "~% nbpoints: ~A" nbpoints)
	       (setf vec3-vector (make-array nbpoints :element-type 'vec3))
	       (incf state))
	      ;; read a  "x y z ... xn yz zn" line and insert parsed vec3s into vec3-vector 
	      ((eq state 4)
	       (if (or (equalp first-token "") (equalp first-token "POLYGONS"))
		   (incf state)
		   (let* ((comp-list (mapcar #'str-coerce-sf (str-tokenize current-line)))
			  (comp-x (sample-list comp-list 2))
			  (comp-y (sample-list (cdr comp-list) 2))
			  (comp-z (sample-list (cddr comp-list) 2))
			  (vec3-lst (mapcar #'vec3-make comp-x comp-y comp-z)))
		     (mapcar (lambda (vec) (setf (aref vec3-vector vec3-vector-mark) vec)
			       (incf vec3-vector-mark))
			     vec3-lst))))
	      ;; parsing over (the rest is not interesting)	      
	      ((eq state 5)
	       (setf parsing-complete t))
	      (t (format t "vtk file parsing failed: state ~A line: ~S" state current-line)
		 (return 0)))))))




;;reads a vtk file represention one 3D object and returns a
;;vector of lists of vertex index contained inside.
;;files must be in unix form no C-m or C-z
;;we assume a parse-vtk file was already performed successfully and do
;;not reiterate standard tests
(defun vtk-parse-polygons (path)
  "return list of polygons from vtk file"
  (let ((get-headtoken-secure (lambda (line) (if (equalp 'eof line) 'eof
						 (str-first-token line))))
	(list-vector nil)
	(parsing-complete nil))
    (with-open-file
	(stream path :direction :input)
      (do* ((current-line (read-line stream nil 'eof)
			  (read-line stream nil 'eof)) ;current line of the stream
	    (first-token (funcall get-headtoken-secure current-line)
			 (funcall get-headtoken-secure current-line))
	    (state 0)
	    (nbpolygons nil)
	    (poly-index 0))
	  ((or (equal current-line 'eof) parsing-complete) list-vector) 

	(cond ((equalp first-token "#") 'comment)
	      ;; read introductory lines up to polygon entry	      	      
	      ((eq state 0)
	       (when (equalp first-token "POLYGONS")
		 (progn 
		   (setf nbpolygons (str-coerce-i (str-fetch-arg "POLYGONS" current-line)))
		   (setf list-vector (make-array nbpolygons :initial-element nil :element-type 'list))
		   (incf state))))
	       
	      ;; read polygons up to first line not starting by integer
	      ((eq state 1)
	       (if (not (typep (read-from-string first-token) 'integer))
		   (incf state)
		   (progn (setf (aref list-vector poly-index)
				(cdr (mapcar #'str-coerce-i (str-tokenize current-line))))
			  (incf poly-index))))
	      
	      ;; parsing over (the rest is not interesting)	      
	      ((eq state 2)
	       (setf parsing-complete t))
	      (t (format t "vtk file parsing failed: state ~A line: ~S" state current-line)
		 (return 0)))))))



(defun vtk-parse-vertices-str (filename)
  "same as vtk-parse-vertices, but file is passed as a string filename"
  (let ((path (make-pathname :name filename)))
    (vtk-parse-vertices path)))


(defun vtk-parse-polygons-str (filename)
  "same as vtk-parse-polygons, but file is passed as a string filename"
  (let ((path (make-pathname :name filename)))
    (vtk-parse-polygons path)))



(defun vtk-parse-tset (file-pattern)
  "returns a list of vec3-vectors representing vtk files matching
given pattern"
  (let* ((str-files (directory file-pattern))
	 (str-lands-vectors (list))) ;;list of vec3-vectors
    (dolist (file str-files (nreverse str-lands-vectors))
      (push (vtk-parse-vertices file) str-lands-vectors))))



(defun write-tset (tset outfile)
  "from a vtk training set (list of vec3-vectors, ie shape instances),
writes a single file featuring all shape vectors (one instance per
line).  returns shape number and shape dimension."
  (flet ((vec3-vector-write
	  (vec3-vector outstream)
	  (progn
	    (map nil 
		 #'(lambda (vec3) (write-string (concatenate 'string (vec3-output-stringize vec3) "  ")
						outstream))
		 vec3-vector)
	    (format outstream "~%"))))
    (let ((path (make-pathname :name outfile)))
      (with-open-file (stream path :direction :output :if-exists :supersede)
	(mapcar #'(lambda (vector) (vec3-vector-write vector stream))
		tset))
      (values (length tset) (length (car tset))))))
	      
  

;;todo: objectname sans extension
(defun vtk-file-to-csv (infile outfile &key (objectname nil))
  "convert infile in vtk format to CSV format and write result to outfile. By default, objectname inferred from filename"
  (let* ((in-path (make-pathname :name infile))
	 (out-path (make-pathname :name outfile))
	 (in-vertices (vtk-parse-vertices in-path))
	 (in-triangles (vtk-parse-polygons in-path))
	 (vertexcount (length in-vertices))
	 (polycount (length in-triangles))
	 (objectname (if objectname objectname infile)))
    (with-open-file (ostream out-path :direction :output :if-exists :supersede)
      (flet ((writestr (string) (write-line string ostream))
	     (vertex-str (num vec3)
			 (mkstr num ", " (vec-x vec3) ", " (vec-y vec3) ", " (vec-z vec3)))
	     (triangle-str (num polylist)
			   (mkstr num ", " (car polylist)  ", " (cadr polylist) ", " (caddr polylist))))
			  
	(progn (writestr "Objects, 1 ")
	       (writestr "Object:, First Vertice:, Vertices:, First Triangle:, Triangles: ")
	       (writestr (mkstr "<" objectname ">, 1, " vertexcount ", 1, " polycount))
	       (writestr "")
	       (writestr (mkstr "Vertices, " vertexcount ))
	       (writestr "Vertice:, X:, Y:, Z: " )
	       (dotimes (i vertexcount 'out)
		 (writestr (vertex-str (1+ i) (aref in-vertices i))))
	       (writestr "")
	       (writestr (mkstr "Triangles, " polycount))
	       (writestr "Triangle:, Side1:, Side2:, Side3: ")
	       (dotimes (i polycount 'out)
		 (writestr (triangle-str (1+ i) (aref in-triangles i))))
	       (writestr "" )
	       'done)))))



(defun vtk-vertex-number (infile)
  "retrieve the number of vertices in given vtk file"
  (vtk-parse-vertices (make-pathname :name infile)
		      :only-nb_vertex t))



;;(bindfun parse-vtkfile vtk-parse-vertices)
;;(bindfun parse-poly-vtkfile vtk-parse-polygons)
;;(bindfun parse-vtkfile-str vtk-parse-vertices-str)
;;(bindfun parse-poly-vtkfile-str vtk-parse-polygons-str)
;;(bindfun parse-vtk-tset vtk-parse-tset)