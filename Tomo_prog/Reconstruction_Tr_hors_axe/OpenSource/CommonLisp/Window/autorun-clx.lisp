;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)
;;
;;This program is free software; you can redistribute it and/or
;;modify it under the terms of the GNU General Public License
;;as published by the Free Software Foundation; either version 2
;;of the License, or (at your option) any later version.
;;
;;This program is distributed in the hope that it will be useful,
;;but WITHOUT ANY WARRANTY; without even the implied warranty of
;;MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;GNU General Public License for more details.
;;
;;You should have received a copy of the GNU General Public License
;;along with this program; if not, write to the Free Software
;;Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.



(load "library:subsystems/clx-library.sparcf")
(use-package :xlib)

;; thanks to cmucl-help
(in-package "XLIB")


(defun get-best-authorization (host display protocol)
  (labels ((read-short (stream &optional (eof-errorp t))
             (let ((high-byte (read-byte stream eof-errorp)))
               (and high-byte
                    (dpb high-byte (byte 8 8) (read-byte stream)))))
           (read-short-length-string (stream)
             (let ((length (read-short stream)))
               (let ((string (make-string length)))
                 (dotimes (k length)
                   (setf (schar string k) (card8->char (read-byte stream))))
                 string)))
           (read-short-length-vector (stream)
             (let ((length (read-short stream)))
               (let ((vector (make-array length :element-type '(unsigned-byte 8))))
                 (dotimes (k length)
                   (setf (aref vector k) (read-byte stream)))
                 vector))))
    ;; Original version didn't handle "localhost" correctly -- SEF.
    (if (string= host "localhost")
        (setq host (machine-instance)))
    (let ((pathname (authority-pathname)))
      (when pathname
        (with-open-file (stream pathname :element-type '(unsigned-byte 8)
                                :if-does-not-exist nil)
          (when stream
            (let* ((host-family (ecase protocol
                                  ((:tcp :internet nil) 0)
                                  ;; The remaining protocols are not really supported -- SEF.
                                  ((:dna :DECnet) 1)
                                  ((:chaos) 2)))
                   (host-address (rest (host-address host host-family))))
              (loop
               (let ((family (read-short stream nil)))
                 (cond ((null family) (return (values "" "")))    ; No useful entry found. -- SEF
                       ((eql family 0)
                        (let* ((address (read-short-length-vector stream))
                               (number (parse-integer (read-short-length-string stream)))
                               (auth-name (read-short-length-string stream))
                               (auth-data (read-short-length-vector stream)))
                          (when (and (= family host-family)
                                     (equal host-address (coerce address 'list))
                                     (= number display)
                                     (string= auth-name "MIT-MAGIC-COOKIE-1"))
                            (return (values auth-name auth-data)))))
                       ;; This is the new case.  The cookie contains a string naming the
                       ;; host, then the display number, auth-name and auth-data. -- SEF
                       ((eql family 256)
                        (let* ((hname (read-short-length-string stream))
                               (number (parse-integer (read-short-length-string stream)))
                               (auth-name (read-short-length-string stream))
                               (auth-data (read-short-length-vector stream)))
                          (when (and (string= hname host)
                                (= number display)
                                (string= auth-name "MIT-MAGIC-COOKIE-1"))
                            (return (values auth-name auth-data)))))))))))))))
