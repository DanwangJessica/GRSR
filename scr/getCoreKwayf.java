import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;


public class getCoreKwayf {
	public static void main(String[] args) throws IOException {
	
		
		FileReader fr1 = new FileReader(args[0]);  //maf file
		BufferedReader br1 = new BufferedReader(fr1);
		
		ArrayList<String> genomeNames=new ArrayList<String>();
		String line=br1.readLine();
		line=br1.readLine();
		String[] aLine=line.split("=");
		String[] sLine;
		String[] temp;
	
		int tot_mult=Integer.parseInt(aLine[3]);
		for(int i=0;i<tot_mult;i++){
			line=br1.readLine();
			temp=line.split("\t\t");
			genomeNames.add(temp[0]);
		}
		//System.out.println(genomeNames.toString()); correct
		fr1.close();
		br1.close();
		
		int istart=1; //this is the start when initializing blocks
		int[] block=new int[tot_mult*4+1];
		
		//Initializing blocks
		block[0]=0; //block id=0
		for(int j=1;j<block.length;j=j+4) {
			block[j]=2; //chr =2 which means this genome doesn't hava this synetic element
			block[j+1]=istart; //start=1
			block[j+2]=1; //size=1
			block[j+3]=1; //strand=1
		}
		istart++;
		//System.out.println(Arrays.toString(block)); correct
		
		int mult=0; //mult of current block
		int score=0; //max length of current block
		int gID=0;
		int count=0;
		FileReader fr = new FileReader(args[0]);  //maf file
		BufferedReader br = new BufferedReader(fr);
		line=br.readLine();
		String[] info;
		int start;
		int end;
     
		while(line!=null){
			line=br.readLine();
			if(line==null) break;
			aLine=line.split("=");
			
			temp=aLine[1].split(" ");
			score=Integer.parseInt(temp[0]);
			mult=Integer.parseInt(aLine[3]);
			
			if(mult!=tot_mult){
				for(int h=0;h<mult;h++) br.readLine();
			}
			else{
				for(int m=0;m<mult;m++){
					line=br.readLine(); //this is the s line
					sLine=line.split("\t\t");
					gID=genomeNames.indexOf(sLine[0]);
					info=sLine[1].split(" ");
					
					if(info[2].contains("+")){
						block[gID*4+1]=1;
						block[gID*4+2]=Integer.parseInt(info[0])+1;
						block[gID*4+3]=Integer.parseInt(info[1]);
						block[gID*4+4]=1;
					}
					else{
						block[gID*4+1]=1;
						end=Integer.parseInt(info[3])-Integer.parseInt(info[0]);
						start=end-Integer.parseInt(info[1])+1;
						block[gID*4+2]=start;
						block[gID*4+3]=Integer.parseInt(info[1]);
						block[gID*4+4]=-1;
					}
				}//end of the for loop for each s_line in a block
				
				System.out.println(Arrays.toString(block).replace(",","").replace("[","").replace("]",""));
				
				//initializing the blocks for next use
				for(int j=1;j<block.length;j=j+4) {
					block[j]=2; //chr =2
					block[j+1]=istart; //start=1
					block[j+2]=1; //size=1
					block[j+3]=1; //strand=1
				}
				istart++;
				count++;
			}//end of the blocks which size is larger than cutoff
			
			line=br.readLine();
		}
		
		//System.out.println("The No. of core blocks is " +count);
		
		
	}

}
