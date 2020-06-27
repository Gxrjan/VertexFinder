using System;

/*
    This class represents a constraint(Inequality)
*/
internal class Inequality
{
    Double[] coefficients;
    Boolean redundant;


    public Inequality(Double[] array)
    {
        this.coefficients = array;
        this.redundant = false;
    }

    public Inequality(Double[] array, bool redundant)
    {
        this.coefficients = array;
        this.redundant = redundant;
    }

    public Inequality(int length)
    {
        this.coefficients = new double[length];
        double[] array = new double[length];
        for (int i = 0; i < length - 1; i++)
        {
            coefficients[i] = 0;
        }
        coefficients[length - 1] = 1;
    }

    public bool is_redundant
    {
        get => this.redundant;

        set => this.redundant = value;
    }

    public int Length
    {
        get => this.coefficients.Length;
    }

    public double this[int key]
    {
        get => this.coefficients[key];
    }


    public override string ToString()
    {
        String result = "";
        foreach (double d in this.coefficients)
        {
            result += d + " ";
        }
        return "[" + result.Trim() + "]";
    }


    public double[] get_coefficients()
    {
        return this.coefficients;
    }
}