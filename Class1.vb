Public Class Sc

    Function MatInv(ma(,)) ' Inverse matrix
        Dim JJ = 0
        Dim J = 0
        Dim L = 0
        Dim K = 0
        Dim F = 0
        Dim D As Double = 0
        Dim EA As Double = 0.0
        Dim i
        i = ma.GetLength(0)
        Dim maa(i, i)
        maa = ma.Clone() ' If I don't clone it will change my inherited matrix ma
        J = 0
        JJ = 0
        Dim MB(i - 1, i - 1) As Double
        While JJ < i
            While J < i
                If J = JJ Then
                    MB(JJ, J) = 1
                Else
                    MB(JJ, J) = 0
                End If
                J = J + 1
            End While
            J = 0
            JJ = JJ + 1
        End While
        JJ = 0
        J = 0

        While JJ < i

            While J < i

                D = 1 / maa(JJ, JJ)
                L = JJ
                While K < i

                    maa(L, K) = maa(L, K) * D
                    MB(L, K) = MB(L, K) * D

                    K = K + 1
                End While

                EA = maa(J, JJ)

                If J <> JJ Then
                    F = 0
                    While F < i
                        maa(J, F) = maa(J, F) - (EA * maa(JJ, F))
                        MB(J, F) = MB(J, F) - (EA * MB(JJ, F))
                        F = F + 1
                    End While

                End If

                J = J + 1

            End While
            F = 0
            EA = 0
            D = 0
            K = 0
            J = 0
            JJ = JJ + 1
        End While

        Return MB

    End Function
    Function MatMult(mc(,), md(,)) 'Matrix multiplication
        Dim r1
        Dim c1
        Dim r2
        Dim c2
        r1 = mc.GetLength(0) - 1
        c1 = mc.GetLength(1) - 1
        r2 = md.GetLength(0) - 1
        c2 = md.GetLength(1) - 1
        Dim c(r1, c2)
        Dim a, b, f As Integer
        Dim g As Double
        a = 0
        b = 0
        f = 0
        g = 0
        While a <= r1
            While b <= c2
                While f <= c1
                    g = g + (mc(a, f) * md(f, b))
                    f = f + 1
                End While
                c(a, b) = g
                b = b + 1
                f = 0
                g = 0

            End While
            g = 0
            b = 0
            f = 0
            a = a + 1
        End While
        Return c
    End Function
    Function MatTrans(ByVal mo(,))
        Dim row As Integer
        Dim col As Integer
        row = mo.GetLength(0) - 1
        col = mo.GetLength(1) - 1
        Dim mmo(col, row)
        Dim aa = 0
        Dim bb = 0
        While bb <= row
            While aa <= col
                mmo(aa, bb) = mo(bb, aa)
                aa = aa + 1
            End While
            aa = 0
            bb = bb + 1
        End While
        Return mmo

    End Function

    Function score(sec(), i)
        'dist is the variable from the optimisation function
        'i is the number of DNAs
        'd is the SCORE oupput array
        'sec() is the actural function output array
        Dim a = 0
        Dim d(i - 1)
        While a < i
            d(a) = (1 / sec(a)) + 0.00001
            a = a + 1

        End While
        a = 0
        Return d
    End Function

    Function MaxScr(d())

        'd is the normalised SCORE Array
        'index is a variable that will find the maximum value in the score array
        'indexx is note used in this function but it would return the array index of the maximum score
        Dim index
        Dim indexx
        index = d.Max()
        indexx = Array.IndexOf(d, index)
        Dim Ace
        Ace = 1 / d.Max()
        Return Ace
    End Function

    Function IndMxScr(d())

        'd is the normalised SCORE Array
        'index is a variable that will find the maximum value in the score array
        'indexx is note used in this function but it would return the array index of the maximum score
        Dim index
        Dim indexx
        index = d.Max()
        indexx = Array.IndexOf(d, index)
        Dim Ace
        Ace = 1 / d.Max()
        Return indexx
    End Function

    Function Opt(a(,) As Double, f() As Double, J As Double, K As Double, pc As Double, Mu As Double, rn1(,) As Double)
        'Notes:
        '1. J is the number of DNAs 
        '2. K is the number of genes in each DNA
        '3. a() is the initial population array
        '4. f() is the initial score array
        '5 pc is crossover rate
        '6 Mu is mutation rate
        '7 rn1 specifies the min/max range of each gene
        '------------------------------------------
        '------------------------------------------0.

        'Caclculate the total of 1/distace {1}
        Dim iii = 0, cc As Double = 0
        While iii < J
            cc = f(iii) + cc
            iii = iii + 1
        End While
        iii = 0
        'END{1}
        ' calculate the probability of each DNA based (1/distance) / total of (1/distance){2}
        Dim nmn(J - 1) As Double
        While iii < J
            nmn(iii) = f(iii) / cc
            iii = iii + 1
        End While
        'END{2}
        'Calculate the cumulative probability from step 2 at each DNA {3}
        iii = 0
        Dim ff(J - 1) As Double
        Dim v As Double = 0
        While iii < J
            ff(iii) = nmn(iii) + v
            v = ff(iii)
            iii = iii + 1
        End While
        'END {3}
        'Step 4 - Create a random array {4}
        iii = 0
        v = 0
        Dim P = 0
        Dim aa = 1
        Dim N(J - 1, K - 1)
        Dim q = 0
        '/RR represents the Random array-----------------------
        '/FF represents the cumulative probability from the population fitness function --------------
        Dim RR(J - 1)
        'Create a random array-----2.3
        Dim ds = 0
        While ds < J
            Randomize()
            RR(ds) = (CInt(Math.Ceiling(Rnd() * 999)) + 0) / 1000
            ds = ds + 1
        End While
        ds = 0
        'END {4}
        'Selection Process {5}
        While P < J
            While aa < J
                If RR(P) <= ff(0) Then
                    While q < K
                        N(P, q) = a(aa - 1, q)
                        q = q + 1

                    End While
                    Exit While
                End If
                If (RR(P) > ff(aa - 1) And RR(P) <= ff(aa)) Then
                    While q < K
                        N(P, q) = a(aa, q)
                        q = q + 1
                    End While

                    aa = J
                    q = 0

                Else

                    aa = aa + 1
                    q = 0

                End If


            End While

            aa = 1
            q = 0
            P = P + 1

        End While
        P = 0
        q = 0
        aa = 1

        'END{5}

        'Selecting Parents for cross-over  {6}
        Dim s = 0
        Dim ss = 0
        Dim M(J - 1, K - 1)
        While s < J
            While ss < K
                M(s, ss) = 998
                ss = ss + 1

            End While
            s = s + 1
        End While
        Dim OP(J - 1)
        Dim Pr = 0

        'OP() is a new random array
        Dim kk = 0
        While kk < J
            Randomize()
            OP(kk) = (CInt(Math.Ceiling(Rnd() * 999)) + 0) / 1000
            If OP(kk) <= pc Then
                While Pr < K
                    M(kk, Pr) = N(kk, Pr)
                    Pr = Pr + 1
                End While
                kk = kk + 1
                Pr = 0
            Else

                ss = 0
                While ss < K
                    M(kk, ss) = 9982016
                    ss = ss + 1

                End While
                kk = kk + 1
                Pr = 0


            End If
        End While
        s = 0
        ss = 0
        kk = 0
        Pr = 0
        'END {6}

        'Cross-Over {7}


        'II is an identifier array that flag the row number of a parent to cross over
        Dim II(J - 1)
        Dim nr = 0
        While nr < J
            If M(nr, 0) <> 9982016 Then
                II(nr) = nr
            Else
                II(nr) = -1

            End If
            nr = nr + 1
        End While
        nr = 0
        'Cross-over operation finalised and transferring the children of parents from matrix M() to Final matrix KO()
        Dim ro = 0
        Dim hh = ro + 1
        Dim co
        Dim fo = 0
        Dim KO(J - 1, K - 1) '/ KO will eventually be replaced with matrix M
        Dim rc(J - 1)
        Dim yyy = 0
        Dim xxx = 0
        While ro < J
            If M(ro, 0) = 9982016 Then

                While yyy < K
                    KO(ro, yyy) = N(ro, yyy)
                    yyy = yyy + 1

                End While

            End If
            xxx = 0
            If ro = J - 1 Then
                yyy = 0

            Else
                yyy = ro + 1
            End If
            If M(ro, 0) <> 9982016 Then
                While yyy < J
                    If M(yyy, 0) = 9982016 Then
                        If (yyy = J - 1) Then
                            yyy = 0
                        Else
                            yyy = yyy + 1


                        End If
                    Else
                        Randomize()
                        co = CInt(Math.Ceiling(Rnd() * (K - 2)))
                        While xxx <= co
                            KO(ro, xxx) = M(ro, xxx)
                            xxx = xxx + 1

                        End While
                        xxx = co + 1
                        While xxx < K
                            KO(ro, xxx) = M(yyy, xxx)
                            xxx = xxx + 1

                        End While
                        xxx = 0
                        yyy = J
                    End If
                End While
            End If
            xxx = 0
            yyy = 0
            ro = ro + 1
        End While

        'END {7}

        'Mutation {8}


        Dim TNGM
        Dim row = KO.GetLength(0)
        Dim Column = KO.GetLength(1)
        If Mu >= 0.1 Then
            TNGM = Math.Round(Mu * (row * Column)) '--- Total number of genetic mutations based on the total number of genes from the population and the mutation rate
        Else
            TNGM = 1
        End If
        Dim YA(TNGM - 1)
        Dim XA(TNGM - 1)
        Dim SA = 0
        Dim r As New Random()


        While SA < TNGM

            ' create 2 matrices that will randomly select the row i and column j of the location of the mutation
            YA(SA) = r.Next(K)
            XA(SA) = r.Next(J)
            SA = SA + 1


        End While
        SA = 0
        Dim rroo As New Random
        'Apply Mutation to the final matrix KO() by changing a genera with a randon number 'r.next()'-----4.1
        While SA < TNGM
            Randomize()
            KO(XA(SA), YA(SA)) = rroo.Next(rn1(0, YA(SA)), rn1(1, YA(SA)))
            SA = SA + 1
        End While
        SA = 0

        'Transfer the final population after mutation to the next generation by transferring the finalised population from KO() to a() -----4.2
        Dim loo = 0, loopl = 0
        While loo < J
            While loopl < K
                a(loo, loopl) = KO(loo, loopl)
                loopl = loopl + 1
            End While
            loopl = 0
            loo = loo + 1
        End While
        Return a

        'ENDD {8}
    End Function
End Class
