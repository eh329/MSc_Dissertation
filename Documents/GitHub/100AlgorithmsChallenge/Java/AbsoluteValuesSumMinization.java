package 100AlgorithmsChallenge;

import java.util.Arrays;
import java.util.Scanner;

public class AbsoluteValuesSum
{
	public static void main(String[] args)
	{
		System.out.println("Enter the size of the array: ");
		Scanner sc = new Scanner(System.in);
		int size = sc.nextInt();
		int[] numbers = new int[size];
		
		for(int i = 0; i < size; i++)
		{
			numbers[i] = sc.nextInt();
		}
		
		Arrays.sort();
		
		if(size % 2 == 0)
		{
			int even_index = (size / 2) - 1;
			System.out.println(numbers[even_index]);
		}
		
		else
		{
			int odd_index = size / 2;
			System.out.println(nmbers[odd_index]);
		}
		
		sc.close();
	}
}