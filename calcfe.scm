#!/usr/bin/env gosh
(use srfi-42)
;;(use file.util)
;;(current-directory)
;;(use gauche.process)
;;(call-with-input-process "ls ../" port->string)
(define-constant temp 298.0d0) ;K


(define-constant cal 4.184) ; J
(define-constant k 1.3806504d-23) ; J/K
(define-constant N 6.02214179d23) ;number
(define (k->R n)
  (define (j->cal n)
    (/. n cal))
  (* N (/. (j->cal n) 1000)))
(define-constant R (k->R k))

(define (au->kcal val)
  (* val 627.509451))

;;Average value is needed to calculate sd. So all values may be read before calclating.
(define (read-file proc file-name)
  (with-input-from-file file-name
    (lambda ()
      (port-for-each
       (lambda (line)
         (proc line))
       read-line))))

(define (foo x)
  (if (<= x -1000)
      (+ x)
      x))

(define (get-energy-from-line line)
  (foo (x->number (list-ref (string-split line #/\s+/) 4))))

(define (get-energy file-name)
  (let ((acc '()))
    (read-file
     (lambda (line) (push! acc (get-energy-from-line line)))
     file-name)
    acc))

(define (average lis)
  (letrec ((rec (lambda (lis tmpsum count) 
                  (if (null? lis)
                      (if (not (= count 0))
                          (/. tmpsum count)
                          (errorf "~s" "zero divide error."))
                      (rec (cdr lis) (+ tmpsum (car lis)) (+ count 1))))))
    (rec lis 0 0)))

(define (var lis)
  (letrec ((rec (lambda (lis av count acc)
                  (if (null? lis)
                      (/ (apply + acc) (- count 1))
                      (rec (cdr lis) av (+ count 1)
                           (cons (let1 val (- (car lis) av) (* val val))
                                 acc))))))
    (rec lis (average lis) 0 '())))

(define (std lis)
  (sqrt (var lis)))

(define (calc-fe ilis jlis)
  ;; dF = -RT ln<exp{=(Ej-Ei)/RT}>
  (letrec ((helper (lambda (i j)
                     (exp (- (/ (au->kcal (- j i)) (* R temp))))))
           (calc-expected-value
            (lambda (i j acc)
              (if (or (null? i) (null? j))
                  (values (average acc) (std acc))
                  (calc-expected-value (cdr i) (cdr j)
                                       (cons
                                        (helper (car i) (car j))
                                        acc))))))
    (receive (av st) (calc-expected-value ilis jlis '())
      ;; ln(x + v) = ln{x(1 + v/x)} = ln(x) + ln(1 + v/x)
      ;; when v/x << 1, ln(1 + v/x) = v/x.
      (define criteria 0.01)
      (let1 v/x (/. st av)
        (if (or (< v/x criteria) (> st av)) ; (- av st) must not be minus.
            (values (* -1 R temp (log av)) v/x (- v/x))
            (let1 center (* -1 R temp (log av))
                    (values
                     center
                     (- (* -1 R temp (log (- av st))) center)
                     (- (* -1 R temp (log (+ av st))) center))))))))

(define (calc-fe-change ifile jfile)
  (calc-fe (get-energy ifile) (get-energy jfile)))

(define (main-routine lis)
  (let ((ifile (ref lis 1))
        (jfile (ref lis 2)))
    (receive (v +s -s)
        (calc-fe-change ifile jfile)
      (format #t "~20s~25s~25s~%" ifile +s -s))
    0))

(define (main args)
  (do-ec (: x 89 141)
         (main-routine (list
                        'dummy
                        (string-append "../result_"
                                       (number->string x)
                                       "-"
                                       (number->string x))
                        (string-append "../result_"
                                       (number->string x)
                                       "-"
                                       (number->string (+ 1 x)))))))


