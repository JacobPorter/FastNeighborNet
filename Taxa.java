package nnet;

import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;
import java.io.Writer;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Vector;
//import jloda.util.Basic;
//import jloda.util.parse.NexusStreamParser;
//import splits.core.SplitsException;
//import splits.core.TaxaSet;

public class Taxa
  extends NexusBlock
{
  private int ntax = 0;
  private Vector taxLabels;
  private Vector taxInfos;
  private boolean mustDetectLabels = false;
  public static final String NAME = "Taxa";
  public static final String FIRSTTAXON = "First-Taxon";
  Taxa originalTaxa;
  TaxaSet hiddenTaxa;
  
  public Taxa()
  {
    clear();
  }
  
  public Taxa(int paramInt)
  {
    clear();
    setNtax(paramInt);
  }
  
  public void clear()
  {
    this.mustDetectLabels = false;
    this.taxLabels = new Vector();
    this.taxInfos = new Vector();
  }
  
  public Object clone()
  {
    Taxa localTaxa = new Taxa();
    try
    {
      localTaxa.setNtax(getNtax());
      for (int i = 1; i <= getNtax(); i++)
      {
        localTaxa.setLabel(i, getLabel(i));
        localTaxa.setInfo(i, getInfo(i));
      }
    }
    catch (Exception localException)
    {
      //Basic.caught(localException);
    }
    return localTaxa;
  }
  
  public boolean equals(Taxa paramTaxa)
  {
    if (getNtax() != paramTaxa.getNtax()) {
      return false;
    }
    for (int i = 1; i <= getNtax(); i++) {
      if (!getLabel(i).equals(paramTaxa.getLabel(i))) {
        return false;
      }
    }
    return true;
  }
  
  public int getNtax()
  {
    return this.ntax;
  }
  
  public void setNtax(int paramInt)
  {
    if (paramInt < 0) {
      paramInt = 0;
    }
    this.ntax = paramInt;
    this.taxLabels.setSize(paramInt + 1);
    this.taxInfos.setSize(paramInt + 1);
  }
  
  public String getLabel(int paramInt)
  {
    return (String)this.taxLabels.get(paramInt);
  }
  
  public void setLabel(int paramInt, String paramString)
    throws SplitsException
  {
    if ((paramInt <= 0) || (paramInt > this.ntax)) {
      throw new SplitsException("index out of range: " + paramInt);
    }
    if (paramString != null)
    {
      String str = ";():\\";
      for (int i = 0; i < str.length(); i++) {
        if (paramString.indexOf("" + str.charAt(i)) != -1) {
          throw new SplitsException("Illegal character '" + str.charAt(i) + "' in taxon label (" + paramInt + "): '" + paramString + "'");
        }
      }
    }
    this.taxLabels.set(paramInt, paramString);
  }
  
  public int indexOf(String paramString)
  {
    if (paramString.equals("First-Taxon")) {
      return 1;
    }
    return this.taxLabels.indexOf(paramString);
  }
  
  public String getInfo(int paramInt)
  {
    return (String)this.taxInfos.get(paramInt);
  }
  
  public void setInfo(int paramInt, String paramString)
    throws SplitsException
  {
    if ((paramInt <= 0) || (paramInt > this.ntax)) {
      throw new SplitsException("index out of range: " + paramInt);
    }
    this.taxInfos.set(paramInt, paramString);
  }
  
  public static void showUsage(PrintStream paramPrintStream)
  {
    paramPrintStream.println("BEGIN TAXA;");
    paramPrintStream.println("DIMENSIONS NTAX=number-of-taxa;");
    paramPrintStream.println("[TAXLABELS taxon_1 taxon_2 ... taxon_ntax;]");
    paramPrintStream.println("[TAXINFO info_1 info_2 ... info_ntax;]");
    paramPrintStream.println("END;");
  }
  
  public void write(Writer paramWriter)
    throws IOException
  {
    paramWriter.write("\nBEGIN Taxa;\n");
    paramWriter.write("DIMENSIONS ntax=" + getNtax() + ";\n");
    paramWriter.write("TAXLABELS\n");
    for (int i = 1; i <= getNtax(); i++) {
      paramWriter.write("[" + i + "] '" + getLabel(i) + "'\n");
    }
    paramWriter.write(";\n");
    if (hasInfos())
    {
      paramWriter.write("TAXINFO\n");
      for (int i = 1; i <= getNtax(); i++) {
        paramWriter.write("[" + i + "] '" + getInfo(i) + "'\n");
      }
      paramWriter.write(";\n");
    }
    paramWriter.write("END; [Taxa]\n");
  }
  
  private boolean hasInfos()
  {
    for (int i = 1; i <= getNtax(); i++) {
      if ((getInfo(i) != null) && (getInfo(i).length() > 0)) {
        return true;
      }
    }
    return false;
  }
  
  public void write(Writer paramWriter, Taxa paramTaxa)
    throws IOException
  {
    write(paramWriter);
  }
  
//  public void read(NexusStreamParser paramNexusStreamParser)
//    throws IOException, SplitsException
//  {
//    paramNexusStreamParser.matchBeginBlock("Taxa");
//    paramNexusStreamParser.matchIgnoreCase("DIMENSIONS ntax=");
//    setNtax(paramNexusStreamParser.getInt());
//    paramNexusStreamParser.matchIgnoreCase(";");
//    int i;
//    if (paramNexusStreamParser.peekMatchIgnoreCase("taxlabels"))
//    {
//      paramNexusStreamParser.matchIgnoreCase("taxlabels");
//      if (paramNexusStreamParser.peekMatchIgnoreCase("_detect_"))
//      {
//        paramNexusStreamParser.matchIgnoreCase("_detect_");
//        setMustDetectLabels(true);
//      }
//      else
//      {
//        for (i = 1; i <= getNtax(); i++)
//        {
//          setLabel(i, paramNexusStreamParser.getLabelRespectCase());
//          for (int j = 1; j < i; j++) {
//            if (getLabel(i).equals(getLabel(j))) {
//              throw new IOException("Line " + paramNexusStreamParser.lineno() + ": " + getLabel(i) + " appears twice in the taxa block");
//            }
//          }
//        }
//      }
//      paramNexusStreamParser.matchIgnoreCase(";");
//    }
//    else
//    {
//      setMustDetectLabels(true);
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("taxinfo"))
//    {
//      paramNexusStreamParser.matchIgnoreCase("taxinfo");
//      for (i = 1; i <= getNtax(); i++) {
//        setInfo(i, paramNexusStreamParser.getLabelRespectCase());
//      }
//      paramNexusStreamParser.matchIgnoreCase(";");
//    }
//    paramNexusStreamParser.matchEndBlock();
//    if (getOriginalTaxa() == null) {
//      setOriginalTaxa((Taxa)clone());
//    }
//  }
//  
  public String toString()
  {
    StringWriter localStringWriter = new StringWriter();
    try
    {
      write(localStringWriter);
    }
    catch (IOException localIOException)
    {
      //Basic.caught(localIOException);
    }
    return localStringWriter.toString();
  }
  
  public TaxaSet getTaxaSet()
  {
    TaxaSet localTaxaSet = new TaxaSet();
    localTaxaSet.set(1, getNtax());
    return localTaxaSet;
  }
  
  public boolean getMustDetectLabels()
  {
    return this.mustDetectLabels;
  }
  
  public void setMustDetectLabels(boolean paramBoolean)
  {
    this.mustDetectLabels = paramBoolean;
    if ((!paramBoolean) && ((getOriginalTaxa() == null) || ((getNtax() > 0) && (getOriginalTaxa().getLabel(1) == null)))) {
      setOriginalTaxa((Taxa)clone());
    }
  }
  
  public static Taxa getInduced(Taxa paramTaxa, TaxaSet paramTaxaSet)
  {
    Taxa localTaxa = new Taxa();
    try
    {
      if (paramTaxaSet != null) {
        localTaxa.setNtax(paramTaxa.getNtax() - paramTaxaSet.cardinality());
      } else {
        localTaxa.setNtax(paramTaxa.getNtax());
      }
      int i = 0;
      for (int j = 1; j <= paramTaxa.getNtax(); j++) {
        if ((paramTaxaSet == null) || (!paramTaxaSet.get(j)))
        {
          i++;
          localTaxa.setLabel(i, paramTaxa.getLabel(j));
          localTaxa.setInfo(i, paramTaxa.getInfo(j));
        }
      }
    }
    catch (SplitsException localSplitsException)
    {
      //Basic.caught(localSplitsException);
    }
    return localTaxa;
  }
  
  public void hideTaxa(TaxaSet paramTaxaSet)
  {
    if ((this.originalTaxa == null) && ((paramTaxaSet == null) || (paramTaxaSet.cardinality() == 0))) {
      return;
    }
    this.hiddenTaxa = paramTaxaSet;
    if (this.originalTaxa == null) {
      this.originalTaxa = ((Taxa)clone());
    }
    Taxa localTaxa = getInduced(this.originalTaxa, paramTaxaSet);
    if (!localTaxa.equals(this)) {
      try
      {
        setNtax(localTaxa.getNtax());
        for (int i = 1; i <= getNtax(); i++)
        {
          setLabel(i, localTaxa.getLabel(i));
          setInfo(i, localTaxa.getInfo(i));
        }
      }
      catch (SplitsException localSplitsException)
      {
        //Basic.caught(localSplitsException);
      }
    }
  }
  
  public Taxa getOriginalTaxa()
  {
    return this.originalTaxa;
  }
  
  public void setOriginalTaxa(Taxa paramTaxa)
  {
    this.originalTaxa = paramTaxa;
  }
  
  public TaxaSet getHiddenTaxa()
  {
    return this.hiddenTaxa;
  }
  
  public void hideAdditionalTaxa(TaxaSet paramTaxaSet)
    throws SplitsException
  {
    if (this.originalTaxa == null) {
      this.originalTaxa = ((Taxa)clone());
    }
    TaxaSet localTaxaSet = new TaxaSet();
    int i = getOriginalTaxa().getNtax();
    if (paramTaxaSet != null) {
      for (int j = 1; j <= i; j++)
      {
        int k = indexOf(this.originalTaxa.getLabel(j));
        if ((k != -1) && (paramTaxaSet.get(k))) {
          localTaxaSet.set(j);
        }
      }
    }
    if ((this.hiddenTaxa != null) && (localTaxaSet != null) && (this.hiddenTaxa.intersects(localTaxaSet))) {
      throw new SplitsException("hidden <" + this.hiddenTaxa + "> and additional <" + localTaxaSet + "> intersect");
    }
    int j = 0;
    for (int k = 1; k <= i; k++) {
      if (((this.hiddenTaxa == null) || (!this.hiddenTaxa.get(k))) && ((paramTaxaSet == null) || (!localTaxaSet.get(k)))) {
        j++;
      }
    }
    setNtax(j);
    j = 0;
    for (int k = 1; k <= i; k++) {
      if (((this.hiddenTaxa == null) || (!this.hiddenTaxa.get(k))) && ((paramTaxaSet == null) || (!localTaxaSet.get(k))) && (getOriginalTaxa() != null))
      {
        setLabel(++j, getOriginalTaxa().getLabel(k));
        setInfo(j, getOriginalTaxa().getInfo(k));
      }
    }
  }
  
  public List getAllLabels()
  {
    LinkedList localLinkedList = new LinkedList();
    for (int i = 1; i <= getNtax(); i++) {
      localLinkedList.add(getLabel(i));
    }
    return localLinkedList;
  }
  
  public static void show(String paramString, Taxa paramTaxa)
  {
    StringWriter localStringWriter = new StringWriter();
    try
    {
      paramTaxa.write(localStringWriter, paramTaxa);
    }
    catch (IOException localIOException)
    {
      //Basic.caught(localIOException);
    }
    System.err.println(paramString + ":\n" + localStringWriter.toString());
  }
  
  public void add(String paramString)
  {
    add(paramString, null);
  }
  
  public void add(String paramString1, String paramString2)
  {
    try
    {
      Taxa localTaxa = (Taxa)clone();
      setNtax(getNtax() + 1);
      for (int i = 1; i <= localTaxa.getNtax(); i++)
      {
        setLabel(i, localTaxa.getLabel(i));
        setInfo(i, localTaxa.getInfo(i));
      }
      setLabel(getNtax(), paramString1);
      setInfo(getNtax(), paramString2);
    }
    catch (SplitsException localSplitsException)
    {
      //Basic.caught(localSplitsException);
    }
  }
  
  public boolean contains(Collection paramCollection)
  {
    Iterator localIterator = paramCollection.iterator();
    while (localIterator.hasNext())
    {
      String str = (String)localIterator.next();
      if (indexOf(str) == -1) {
        return false;
      }
    }
    return true;
  }
  
  public Set<String> getLabels(TaxaSet paramTaxaSet)
  {
    HashSet localHashSet = new HashSet();
    for (int i = paramTaxaSet.getBits().nextSetBit(0); i != -1; i = paramTaxaSet.getBits().nextSetBit(i + 1)) {
      localHashSet.add(getLabel(i));
    }
    return localHashSet;
  }
}
