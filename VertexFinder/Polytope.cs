using System;
using System.Text;
using System.Collections;
using System.Collections.Generic;
using Google.OrTools.LinearSolver;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Threading;
using MathNet.Numerics;

/*
    This class represents a Polytope
*/
internal class Polytope
{
    // Needed for Multithreading
    private readonly object _sync = new object();

    // Represents the constaints of a polytope
    public Inequality[] constraints;

    //This variable is the way we defferentiate two polytopes.
    //Plus it's the way we see if i'th inequality is transformed to equality
    public string id_string;

    // Array of facets. Direct children of this polytope in a face lattice
    public Polytope[] facets;

    // 0 <= 1
    Inequality shadow;

    // This variable contains all faces of the "original"(the one read from the file) polytope.
    // This variable is "shared" among all faces of a face lattice. Meaning that both original polytope and it's vertices
    // see the same thing
    public HashSet<Polytope> faces;


    // This variable represents the number of times explore function has been called on this polytope
    private int visitCount;


    // This variable represents the number of iteration, that it has been discovered
    public int discoveryIteration = -1;


    // A set of vertices of the original polytope.
    public HashSet<Polytope> vertices = new HashSet<Polytope>();

    public Polytope(Inequality[] constraints, string parent_id_string, HashSet<Polytope> exploredFaces)
    {
        this.constraints = constraints;
        this.id_string = parent_id_string;
        this.id_string = calculate_string(this.constraints, parent_id_string);
        this.shadow = new Inequality(constraints[0].Length);
        this.faces = exploredFaces;
        this.visitCount = 0;
    }


    public Polytope(Inequality[] constraints)
    {
        this.constraints = constraints;
        this.id_string = new string('0', constraints.Length);
        this.shadow = new Inequality(constraints[0].Length);
        faces = new HashSet<Polytope>();
        faces.Add(this);
        this.visitCount = 0;
    }

    public int visit_count
    {
        get => this.visitCount;
        set => this.visitCount = value;
    }


    /// <summary>
    /// This method considres an id_string of a polytope and figures out which inequalities should be equalities.
    /// </summary>
    /// <param name="constraints">Ineqaultites of a polytope</param>
    /// <param name="old_id_string">Previsous id_string</param>
    /// <returns></returns>
    private string calculate_string(Inequality[] constraints, string old_id_string)
    {
        StringBuilder builder = new StringBuilder(old_id_string);
        for (int i = 0; i < constraints.Length; i++)
        {
            if (old_id_string[i] == '0')
            {
                if (should_be_equality(i))
                {
                    builder[i] = '1';
                }
            }
        }
        return builder.ToString();
    }



    /// <summary>
    /// Determines if i'th inequality is redundant or not
    /// </summary>
    /// <param name="inequality_index">The index of inequality</param>
    /// <returns></returns>
    public bool is_redundant(int inequality_index)
    {
        if (this.constraints[inequality_index].is_redundant)
            return true;

        Inequality[] inequalities = get_non_redundant_inequalitites_except(inequality_index);
        if (inequalities.Length == 1) // means that there is only shadow inequality in there
            return false;

        Inequality suspect = this.constraints[inequality_index];


        Solver solver = Solver.CreateSolver("redundantProgram", "GLOP_LINEAR_PROGRAMMING");
        Variable[] variables = new Variable[inequalities.Length]; // create variable for each vector
        for (int i = 0; i < inequalities.Length; i++)
        {
            variables[i] = solver.MakeNumVar(0, double.PositiveInfinity, "var " + i); //make them non-negative(Conical combination)
        }

        for (int i = 0; i < suspect.Length; i++)
        {
            Constraint ct = solver.MakeConstraint(suspect[i], suspect[i], "ct " + i);
            for (int j = 0; j < variables.Length; j++)
            {
                ct.SetCoefficient(variables[j], inequalities[j][i]);
            }
        }
        Solver.ResultStatus result = solver.Solve();
        bool b = result == Solver.ResultStatus.FEASIBLE || result == Solver.ResultStatus.UNBOUNDED ||
               result == Solver.ResultStatus.OPTIMAL;
        this.constraints[inequality_index].is_redundant = b;
        return b;
    }



    /// <summary>
    /// Returns non redundant inequalities except i'th
    /// </summary>
    /// <param name="inequality_index">the index of inuequality we shouldn't include</param>
    /// <returns></returns>
    private Inequality[] get_non_redundant_inequalitites_except(int inequality_index)
    {
        List<Inequality> result = new List<Inequality>();
        for (int i = 0; i < this.constraints.Length; i++)
        {
            if (i == inequality_index)
                continue;
            if (id_string[i] == '0')
            {
                if (!constraints[i].is_redundant)
                {
                    result.Add(constraints[i]);
                }
            }
            else
            {
                result.Add(this.constraints[i]);
                result.Add(flip(this.constraints[i]));
            }
        }
        result.Add(shadow);
        return result.ToArray();
    }



    /// <summary>
    /// Takes an inequality and returns the same inequalities with coefficients multiplied by -1
    /// </summary>
    /// <param name="inequality">Inequality to be flipped</param>
    /// <returns></returns>
    private Inequality flip(Inequality inequality)
    {
        double[] newCoefficients = new double[inequality.Length];
        double[] oldCoefficients = inequality.get_coefficients();
        for (int i = 0; i < newCoefficients.Length; i++)
        {
            newCoefficients[i] = oldCoefficients[i] * -1;
        }
        return new Inequality(newCoefficients);
    }



    /// <summary>
    /// Returns an array with all inequalitites
    /// </summary>
    /// <returns>Returns an array with all inequalitites</returns>
    private Inequality[] get_all_inequalities()
    {
        List<Inequality> result = new List<Inequality>();
        for (int i = 0; i < this.constraints.Length; i++)
        {
            if (id_string[i] == '0')
            {
                result.Add(this.constraints[i]);
            }
        }
        return result.ToArray();
    }



    /// <summary>
    /// Returns an array of all equalitites except i'th
    /// </summary>
    /// <param name="equality_index">the index of inuequality we shouldn't include</param>
    /// <returns></returns>
    private Inequality[] get_equalities_except(int equality_index)
    {
        List<Inequality> result = new List<Inequality>();
        for (int i = 0; i < this.constraints.Length; i++)
        {
            if (id_string[i] == '1' && i != equality_index)
            {
                result.Add(this.constraints[i]);
            }
        }
        return result.ToArray();
    }



    /// <summary>
    /// Determines if i'th equality should be equality
    /// </summary>
    /// <param name="inequality_index">The index of an ineqaulity</param>
    /// <returns>True if should be equality</returns>
    public bool should_be_equality(int inequality_index)
    {
        Inequality[] equalities = get_equalities_except(inequality_index);
        if (equalities.Length == 0)
            return false;
        Inequality suspect = this.constraints[inequality_index];
        double[][] columnArray = new double[equalities.Length][];
        double[][] columnArray2 = new double[equalities.Length + 1][];
        for (int i = 0; i < columnArray.Length; i++)
        {
            columnArray[i] = equalities[i].get_coefficients();
            columnArray2[i] = equalities[i].get_coefficients();
        }
        columnArray2[columnArray.Length] = suspect.get_coefficients();



        Matrix<double> A = DenseMatrix.OfColumnArrays(columnArray);
        Matrix<double> B = DenseMatrix.OfColumnArrays(columnArray2);
        bool result = A.Rank() == B.Rank();
        return result;
    }



    /// <summary>
    /// Computes the facets of a polytope
    /// </summary>
    /// <returns>Array of polytopes(facets)</returns>
    public Polytope[] get_facets()
    {
        List<Polytope> result = new List<Polytope>();
        for (int i = 0; i < constraints.Length; i++)
        {
            if (id_string[i] != '1' && !constraints[i].is_redundant)
            {
                if (is_redundant(i))
                {
                    constraints[i].is_redundant = true;
                }
                else
                {
                    StringBuilder builder = new StringBuilder(id_string);
                    builder[i] = '1';
                    string new_id_string = builder.ToString();
                    Polytope facet = new Polytope(clone_constraints(), new_id_string, faces);
                    Polytope newFacet;
                    if (!this.faces.Contains(facet))
                    {
                        this.faces.Add(facet);
                    }
                    else
                    {
                        this.faces.TryGetValue(facet, out facet);
                    }

                    result.Add(facet);
                }
            }
        }
        return result.ToArray();
    }


    /// <summary>
    /// Creates a new pointer to inequalites and clones all inequalites from this polytope and puts it under the new pointer
    /// </summary>
    /// <returns>Array of inequalites</returns>
    public Inequality[] clone_constraints()
    {
        List<Inequality> result = new List<Inequality>();
        foreach (Inequality ineq in constraints)
        {
            Inequality newIneq = new Inequality(ineq.get_coefficients(), ineq.is_redundant);
            result.Add(newIneq);
        }

        return result.ToArray();
    }


    public override string ToString()
    {
        string result = "";
        foreach (Inequality ineq in constraints)
        {
            result += ineq + " Is redundant?: " + ineq.is_redundant + "\n";
        }
        result += "id_string: " + id_string + "\n";
        result += "Is vertex?: " + this.is_vertex();
        return result;
    }



    /// <summary>
    /// Determines if this is a vertex polytope
    /// </summary>
    /// <returns>True if it is a vertex</returns>
    public bool is_vertex()
    {
        for (int i = 0; i < constraints.Length; i++)
        {
            if (id_string[i] == '0' && !is_redundant(i))
            {
                return false;
            }
        }
        return true;
    }


    /// <summary>
    /// Returns all inequalitites that are equalitites
    /// </summary>
    /// <returns>Returns all inequalitites that are equalitites</returns>
    private Inequality[] get_all_equalities()
    {
        List<Inequality> result = new List<Inequality>();
        for (int i = 0; i < this.constraints.Length; i++)
        {
            if (id_string[i] == '1')
            {
                result.Add(this.constraints[i]);
            }
        }
        return result.ToArray();
    }


    /// <summary>
    /// Performs a simple dfs on a face lattice
    /// </summary>
    /// <param name="visited">List of polytopes that dfs already visited</param>
    public void dfs(List<Polytope> visited)
    {
        visited.Add(this);
        if (is_vertex())
        {

            calculate_point_and_print();
        }

        if (this.facets == null)
            this.facets = get_facets();
        foreach (Polytope facet in facets)
            if (!visited.Contains(facet))
                facet.dfs(visited);

    }



    /// <summary>
    /// When some vertex is reached we can caluclate the vector corresponding to this vertex.
    /// After calculating we display it in a format: "{discovery iteration}:{vector}"
    /// </summary>
    private void calculate_point_and_print()
    {
        Inequality[] equalities = get_all_equalities();
        double[][] matrixCoeff = new double[equalities.Length][];
        double[] bCoeff = new double[equalities.Length];
        for (int i = 0; i < equalities.Length; i++)
        {
            matrixCoeff[i] = new double[equalities[i].Length - 1];
            for (int j = 0; j < equalities[i].Length - 1; j++)
            {
                matrixCoeff[i][j] = equalities[i][j];
            }
            bCoeff[i] = equalities[i][equalities[i].Length - 1];
        }
        Matrix<double> A = DenseMatrix.OfRowArrays(matrixCoeff);
        Vector<double> b = DenseVector.OfArray(bCoeff);
        Vector<double> x = A.Solve(b);
        Console.WriteLine(x);

    }



    /// <summary>
    /// If this polytope is a vertex that we can calculate the coordinates of it.
    /// </summary>
    /// <returns></returns>
    public Vector<double> calculate_point()
    {
        Inequality[] equalities = get_all_equalities();
        double[][] matrixCoeff = new double[equalities.Length][];
        double[] bCoeff = new double[equalities.Length];
        for (int i = 0; i < equalities.Length; i++)
        {
            matrixCoeff[i] = new double[equalities[i].Length - 1];
            for (int j = 0; j < equalities[i].Length - 1; j++)
            {
                matrixCoeff[i][j] = equalities[i][j];
            }
            bCoeff[i] = equalities[i][equalities[i].Length - 1];
        }
        Matrix<double> A = DenseMatrix.OfRowArrays(matrixCoeff);
        Vector<double> b = DenseVector.OfArray(bCoeff);
        Vector<double> x = A.Solve(b);
        return x;
    }


    public static bool operator ==(Polytope p1, Polytope p2)
    {
        return p1.id_string == p2.id_string;
    }


    public static bool operator !=(Polytope p1, Polytope p2)
    {
        return p1.id_string != p2.id_string;
    }


    public override int GetHashCode()
    {
        return this.id_string.GetHashCode();
    }


    public override bool Equals(object obj)
    {
        return this == (Polytope)obj;
    }



    /// <summary>
    /// This recursive function explores a facce lattice by picking each facet randomly or by picking the least visited facet.
    /// </summary>
    /// <param name="iteration">The number of times "explore" has been called on the original polytope</param>
    /// <param name="vertices">Set of vertices that we have already dicovered</param>
    public void explore(int iteration, HashSet<Polytope> vertices)
    {
        this.visitCount++;
        if (this.discoveryIteration == -1)
            this.discoveryIteration = iteration;
        if (is_vertex())
        {
            if (!vertices.Contains(this))
            {
                vertices.Add(this);
                print_vertex(vertices.Count);
            }
            return;
        }
        Polytope facet = get_least_visited_facet();
        //Polytope facet = get_random_facet();
        facet.explore(iteration, vertices);

    }

    /// <summary>
    /// Print vertex info in the format "{discovery iteration}:{vector}"
    /// </summary>
    public void print_vertex(int num)
    {
        Vector<double> x = calculate_point();
        Console.Write(num + " " + this.discoveryIteration + ":");
        foreach (double d in x)
        {
            double digit = d;
            if (d.AlmostEqual(0))
                digit = Math.Abs(d);
            Console.Write(String.Format("{0:0.##} ", digit));
        }
        Console.WriteLine();
    }


    /// <summary>
    /// Returns random facet of a polytope
    /// </summary>
    /// <returns></returns>
    public Polytope get_random_facet()
    {
        if (this.facets == null)
            this.facets = get_facets();
        Random rand = new Random();
        return this.facets[rand.Next(this.facets.Length)];
    }



    /// <summary>
    /// Returns the least visited facet of a polytope
    /// </summary>
    /// <returns></returns>
    public Polytope get_least_visited_facet()
    {
        if (this.facets == null)
            this.facets = get_facets();
        Polytope result = facets[0];
        foreach (Polytope facet in facets)
            if (facet.visit_count <= result.visit_count)
                result = facet;
        return result;
    }


    /// <summary>
    /// Outer function of dfs
    /// </summary>
    public void explore_all_faces()
    {
        List<Polytope> visited = new List<Polytope>();
        this.dfs(visited);
    }


    /// <summary>
    /// Every 4 seconds displays info about vertices that have been discovered so far.
    /// Not used.
    /// </summary>
    public void display_vertices()
    {
        while (true)
        {
            lock (_sync)
            {
                Console.Clear();
                Console.WriteLine("Vertices total: " + this.vertices.Count);
                foreach (Polytope vertex in this.vertices)
                {
                    Console.WriteLine(/*vertex.id_string +*/ "        visit count: " + vertex.visit_count);
                    Vector<Double> point = vertex.calculate_point();
                    Console.WriteLine("discovery iteration: " + vertex.discoveryIteration);
                    Console.Write("             Vector: ");
                    foreach (Double d in point)
                    {
                        double digit = d;
                        if (d == 0 || d == -0)
                            digit = Math.Abs(d);
                        Console.Write(String.Format("{0:0.##} ", d));
                    }
                    Console.WriteLine();
                    Console.WriteLine();
                }
            }
            Thread.Sleep(4000);
        }
    }
}
