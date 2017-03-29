
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class getMaxL {

	public static void main(String[] args) throws IOException {
		//args[0] is the input position file start end which can be empty
		//print the max length of the input file
		FileReader fr = new FileReader(args[0]); 
		BufferedReader br = new BufferedReader(fr);
		String line;
		String[] tokens;
		int maxL=0;
		int curL=0;
		while(true){
			line=br.readLine();
			if(line==null) break;
			tokens=line.split(" ");
			curL=Integer.parseInt(tokens[1])-Integer.parseInt(tokens[0])+1;
			if(curL>maxL)
				maxL=curL;
		}
		System.out.println(maxL);
		
		br.close();

	}

}
