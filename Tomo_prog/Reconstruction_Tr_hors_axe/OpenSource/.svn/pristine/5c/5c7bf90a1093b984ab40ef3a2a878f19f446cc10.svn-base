;;; ----------------------------------------------------------------------------
;;; $RCSfile: IOStream.lisp,v $
;;; $Author: hubert $
;;; $Date: 2000/05/05 06:13:58 $
;;; $Revision: 1.1.1.1 $
;;; ----------------------------------------------------------------------------

;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;;(load "/users/these/bailleul/These/OpenSource/CommonLisp/Tools.sparcf")


;;;-----------------------------------------------------------------------------
;;; Tool functions from Paul Graham
;;;-----------------------------------------------------------------------------


;;Paul Graham's 'Ansi Common Lisp'
(defun tokens (str test start)
  (let ((p1 (position-if test str :start start)))
    (if p1
	(let ((p2 (position-if #'(lambda (c) 
				   (not (funcall test c)))
			       str :start p1)))
	  (cons (subseq str p1 p2)
		(if p2 
		    (tokens str test p2) 
		    nil)))
	nil)))


;;;-----------------------------------------------------------------------------
;;; String Functions
;;;-----------------------------------------------------------------------------


(defun str-exclude (x)
  "return nil if the argument is a character which should not appear in the
final string."
  (or (eq #\SPACE x)
      (eq #\" x)
      (eq #\, x)))


(defun str-tokenize (string)
  "return a list of tokens, ie #\SPACE-separated elements, from given
string.  Notice that encountered \"strings\" are included without
double-quotes."
  (assert (not (eq string 'eof)))
  (let ((nospc (lambda (x) (not (str-exclude x)))))
    (tokens string nospc 0)))


(defun str-first-token (string)
  "return the first encountered token in the given string"
  (string-right-trim
   '(#\SPACE)
   (block run
     (do* ((i 0 (1+ i))
	   (c nil)
	   (start nil)
	   (len (length string))
	   (result (make-array len :element-type 'character
		 	       :initial-element #\ ))
	   (result-pos 0)) 
	 ((>= i len) result)
       (setf c (char string i))
       (cond ((str-exclude c) (if start (return-from run result)))
	     (t (setf start t)
		(setf (char result result-pos) c)
		(incf result-pos)))))))


(defun str-fetch-arg (com string)
  "return from the given string the token just after the given one"
  (let ((com-len (length com))
	(com-pos (search com string)))
    (cond (com-pos (str-first-token (subseq string (+ com-len com-pos))))
	  (t nil))))
      
      
(defmacro str-append-to (str app)
  "appends the app string to the given string (side-effect)" 
  `(setf ,str (concatenate 'string ,str ,app)))


(defun str-append (str app)
  "returns the concatenation of given strings"
  (concatenate 'string str app))


(defun str-coerce-sf (str)
  "returns a single-float from given string if possible"
  (coerce (read-from-string str)
	  'single-float))


(defun str-coerce-r (str)
  "returns a real from given string if possible"
  (coerce (read-from-string str)
	  'real))


(defun str-coerce-i (str)
  "returns an integer from given string if possible (if numeric)"
  (coerce (round (read-from-string str))
	  'integer))


(defun l-coerce-sf (values-lst)
  "coerce values of given list to single-float"
  (mapcar #'(lambda (x) (coerce x 'single-float)) values-lst))


(defun l-coerce-df (values-lst)
  "coerce values of given list to double-float"
  (mapcar #'(lambda (x) (coerce x 'double-float)) values-lst))


;;;-----------------------------------------------------------------------------
;;; Stream seeking Functions
;;;-----------------------------------------------------------------------------


(defun stream-seek (pattern stream current-line) 
  "Seeks in stream a given string pattern, if not present in the current-line string.
Returns the line of the stream where the pattern first occurs, NIL if end-of-line reached. 
Side-effet on the given stream."
  (do ((line current-line (read-line stream nil 'eof)))
      ((eq line 'eof) nil)
    (if (search pattern line) (return line))))


(defun stream-multiline-seek (start-pattern stop-pattern-list stream current-line)
  "Extension to seek-in-stream: seeks between a start pattern and any stop-pattern,
commencing from the given current-line. Returns the new current-line as a first value,
the sought string as the second value."
  (do* ((line (stream-seek start-pattern stream current-line)
	      (read-line stream nil 'eof))
	(acc nil)
	(str-seek (lambda (pat) (search pat line))))
      ((eq line 'eof) (values nil acc))
    (if (some #'identity (mapcar str-seek stop-pattern-list))
	(return (values line acc)))
    (str-append-to acc line)))


;;;-----------------------------------------------------------------------------
;;; File utilities
;;;-----------------------------------------------------------------------------


(defun file-to-token-lists (input-file-name)
  "reads a test file: return a list of lines, ie a sublist of string tokens"
  (with-open-file (stream (make-pathname :name input-file-name)
			  :direction :input)
    (loop for line = (read-line stream nil)
	  while line
	  collect (str-tokenize line))))


(defun lists-to-strings (listoflists &key (token-separator " "))
  "from a list of sublists -containing anything- , infer a list of strings"
  (mapcar #'(lambda (lst) (apply #'concatenate 'string
				 (mapcar #'(lambda (elem) (mkstr elem token-separator))
					 lst)))
	  listoflists))



(defmacro list-to-string (lst &key (token-separator " "))
  "converts a list to a string. built upon lists-to-strings"
  `(lists-to-strings (list ,lst) :token-separator ,token-separator))


(defmacro lists-to-string (listsoflists &key (token-separator " "))
  "converts several lists to a single string. built upon lists-to-strings"
  `(apply #'concatenate 'string
	  (lists-to-strings ,listsoflists :token-separator ,token-separator)))


;; !!to use in most cases
(defun stringlist-to-file (string-list output-file-name)
  "writes a list of strings into a file"
  (with-open-file (ostream (make-pathname :name output-file-name)
			   :direction :output :if-exists :supersede)
    (flet ((writestr (string) (write-line string ostream)))
      (mapcar #'writestr string-list))))


(defun string-to-file (in-string output-file-name)
  "writes a string into a file"
  (with-open-file (ostream (make-pathname :name output-file-name)
			   :direction :output :if-exists :supersede)
    (write-line in-string ostream)))
    
    
(defun list-to-file (mylist output-file-name)
  "writes each list element onto a single file line"
  (with-open-file (ostream (make-pathname :name output-file-name)
			   :direction :output :if-exists :supersede)
    (flet ((write-obj (any) (write-line (format nil "~A" any) ostream)))
      (mapcar #'write-obj mylist))))

(defun parse-filename (in-name line-process)
  "reads filename-designated file line per line: applies line-process
function to each string line."
  (with-open-file (stream (make-pathname :name in-name) :direction :input)
    (loop for current-line = (read-line stream nil)
	  while current-line
	  do (funcall line-process current-line))))
