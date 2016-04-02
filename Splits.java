package nnet;

import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.lang.reflect.Array;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import nnet.SplitsSet;
import nnet.TaxaSet;
//import jloda.util.Alert;
//import jloda.util.Basic;
//import jloda.util.Interval;
//import jloda.util.parse.NexusStreamParser;
//import splits.core.SplitsException;
//import splits.core.SplitsSet;
//import splits.core.TaxaSet;

public class Splits
  extends NexusBlock
  implements Cloneable
{
  public static final String NAME = "Splits";
  private Format format = null;
  private Properties properties = null;
  private float threshold = 0.0F;
  private int[] cycle = null;
  private SplitsSet splits = new SplitsSet();
  private Taxa previousTaxa;
  private Splits originalSplits;
  
  public SplitsSet getSplitsSet()
  {
    return this.splits;
  }
  
  public void addSplitsSet(SplitsSet paramSplitsSet)
  {
    for (int i = 1; i <= paramSplitsSet.getNsplits(); i++) {
      add(paramSplitsSet.getSplit(i), paramSplitsSet.getWeight(i));
    }
  }
  
  public Splits() {}
  
  public Splits(int paramInt)
  {
    this();
    this.splits.setNtax(paramInt);
  }
  
  public Format getFormat()
  {
    return this.format;
  }
  
  public Properties getProperties()
  {
    return this.properties;
  }
  
  public float getThreshold()
  {
    return this.threshold;
  }
  
  public void setThreshold(float paramFloat)
  {
    this.threshold = paramFloat;
  }
  
  public int getNtax()
  {
    return this.splits.getNtax();
  }
  
  public void setNtax(int paramInt)
  {
    this.splits.setNtax(paramInt);
  }
  
  public int getNsplits()
  {
    return this.splits.getNsplits();
  }
  
  public void setNsplits(int paramInt)
  {
    this.splits.setNsplits(paramInt);
  }
  
  public void setCycle(int[] paramArrayOfInt)
    throws SplitsException
  {
    if (paramArrayOfInt != null)
    {
      BitSet localBitSet = new BitSet();
      for (int i = 1; i < paramArrayOfInt.length; i++)
      {
        if (localBitSet.get(paramArrayOfInt[i])) {
          throw new SplitsException("setCycle(): Multiple occurence of taxon " + i);
        }
        localBitSet.set(paramArrayOfInt[i]);
      }
    }
    this.cycle = paramArrayOfInt;
  }
  
  public int[] getCycle()
  {
    return this.cycle;
  }
  
  public String getLabel(int paramInt)
  {
    return this.splits.getLabel(paramInt);
  }
  
  public void setLabel(int paramInt, String paramString)
  {
    this.splits.setLabel(paramInt, paramString);
  }
  
  public void clear()
  {
    this.splits.setNtax(0);
    this.splits.setNsplits(0);
    this.splits.setLabels(new Vector());
    this.splits.setSplits(new Vector());
    this.splits.setWeights(new Vector());
  }
  
  public void add(TaxaSet paramTaxaSet)
  {
    this.splits.add(paramTaxaSet);
  }
  
  public void add(TaxaSet paramTaxaSet, float paramFloat)
  {
    this.splits.add(paramTaxaSet, paramFloat);
  }
  
  public void add(TaxaSet paramTaxaSet, float paramFloat1, float paramFloat2)
  {
    this.splits.add(paramTaxaSet, paramFloat1, paramFloat2);
  }
  
  public void add(TaxaSet paramTaxaSet, float paramFloat, String paramString)
  {
    this.splits.add(paramTaxaSet, paramFloat, paramString);
  }
  
  public void add(TaxaSet paramTaxaSet, float paramFloat1, float paramFloat2, String paramString)
  {
    this.splits.add(paramTaxaSet, paramFloat1, paramFloat2, paramString);
  }
  
  public void add(TaxaSet paramTaxaSet, float paramFloat1, float paramFloat2, Interval paramInterval, String paramString)
  {
    this.splits.add(paramTaxaSet, paramFloat1, paramFloat2, paramInterval, paramString);
  }
  
  public void remove(int paramInt)
  {
    this.splits.remove(paramInt);
  }
  
  public TaxaSet get(int paramInt)
  {
    return this.splits.getSplit(paramInt);
  }
  
  public float getWeight(int paramInt)
  {
    return this.splits.getWeight(paramInt);
  }
  
  public float getConfidence(int paramInt)
  {
    return this.splits.getConfidence(paramInt);
  }
  
  public void setConfidence(int paramInt, float paramFloat)
  {
    this.splits.setConfidence(paramInt, paramFloat);
  }
  
  public Interval getInterval(int paramInt)
  {
    return this.splits.getInterval(paramInt);
  }
  
  public void setInterval(int paramInt, Interval paramInterval)
  {
    this.splits.setInterval(paramInt, paramInterval);
  }
  
  public String toLogString(int paramInt)
  {
    return " [" + getLabel(paramInt) + "][" + get(paramInt).toString() + "][" + getWeight(paramInt) + "]";
  }
  
  public String toLogString()
  {
    StringWriter localStringWriter = new StringWriter();
    try
    {
      write(localStringWriter, getNtax());
    }
    catch (IOException localIOException)
    {
      localIOException.printStackTrace();
    }
    return toString() + "\n" + localStringWriter.toString();
  }
  
  public void setWeight(int paramInt, float paramFloat)
  {
    this.splits.setWeight(paramInt, paramFloat);
  }
  
  public int indexOf(String paramString)
  {
    return this.splits.indexOf(paramString);
  }
  
  public static void showUsage(PrintStream paramPrintStream)
  {
    paramPrintStream.println("BEGIN Splits;");
    paramPrintStream.println("\t[DIMENSIONS [NTAX=number-of-taxa] [NSPLITS=number-of-splits];]");
    paramPrintStream.println("\t[FORMAT");
    paramPrintStream.println("\t    [LABELS={LEFT|NO}]");
    paramPrintStream.println("\t    [WEIGHTS={YES|NO}]");
    paramPrintStream.println("\t    [CONFIDENCES={YES|NO}]");
    paramPrintStream.println("\t    [INTERVALS={YES|NO}]");
    paramPrintStream.println("\t;]");
    paramPrintStream.println("\t[THRESHOLD=non-negative-number;]");
    paramPrintStream.println("\t[PROPERTIES");
    paramPrintStream.println("\t    [FIT=non-negative-number]");
    paramPrintStream.println("\t    [leastsquares]");
    paramPrintStream.println("\t    [{COMPATIBLE|CYCLIC|WEAKLY COMPATIBLE|INCOMPATIBLE]");
    paramPrintStream.println("\t;]");
    paramPrintStream.println("\t[CYCLE [taxon_i_1 taxon_i_2 ... taxon_i_ntax];]");
    paramPrintStream.println("\t[SPLITSLABELS label_1 label_2 ... label_nsplits;]");
    paramPrintStream.println("\tMATRIX");
    paramPrintStream.println("\t    [label_1] [weight_1] [confidence_1] split_1,");
    paramPrintStream.println("\t    [label_2] [weight_2] [confidence_2] split_2,");
    paramPrintStream.println("\t    ....");
    paramPrintStream.println("\t    [label_nsplits] [weight_nsplits] [confidence_nsplits] split_nsplits,");
    paramPrintStream.println("\t;");
    paramPrintStream.println("END;");
  }
  
  public void write(Writer paramWriter, Taxa paramTaxa)
    throws IOException
  {
    write(paramWriter, paramTaxa.getNtax());
  }
  
  public void write(Writer paramWriter, int paramInt)
    throws IOException
  {
    paramWriter.write("\nBEGIN Splits;\n");
    paramWriter.write("DIMENSIONS ntax=" + getNtax() + " nsplits=" + getNsplits() + ";\n");
    paramWriter.write("FORMAT");
    if (getFormat().getLabels() == true) {
      paramWriter.write(" labels=left");
    } else {
      paramWriter.write(" labels=no");
    }
    if (getFormat().getWeights() == true) {
      paramWriter.write(" weights=yes");
    } else {
      paramWriter.write(" weights=no");
    }
    if (getFormat().getConfidences() == true) {
      paramWriter.write(" confidences=yes");
    } else {
      paramWriter.write(" confidences=no");
    }
    if (getFormat().getIntervals() == true) {
      paramWriter.write(" intervals=yes");
    } else {
      paramWriter.write(" intervals=no");
    }
    paramWriter.write(";\n");
    if (getThreshold() != 0.0F) {
      paramWriter.write("THRESHOLD=" + getThreshold() + "; \n");
    }
    paramWriter.write("PROPERTIES fit=" + getProperties().getFit());
    if (getProperties().isLeastSquares()) {
      paramWriter.write(" leastsquares");
    }
    switch (getProperties().getCompatibility())
    {
    case 1: 
      paramWriter.write(" compatible");
      break;
    case 2: 
      paramWriter.write(" cyclic");
      break;
    case 3: 
      paramWriter.write(" weakly compatible");
      break;
    case 4: 
      paramWriter.write(" non compatible");
    }
    paramWriter.write(";\n");
    int j;
    if (getCycle() != null)
    {
      paramWriter.write("CYCLE");
      int[] arrayOfInt = getCycle();
      for (j = 1; j < Array.getLength(arrayOfInt); j++) {
        paramWriter.write(" " + arrayOfInt[j]);
      }
      paramWriter.write(";\n");
    }
    paramWriter.write("MATRIX\n");
    for (int i = 1; i <= getNsplits(); i++)
    {
      j = Math.min(get(i).cardinality(), paramInt - get(i).cardinality());
      paramWriter.write("[" + i + ", size=" + j + "]" + " \t");
      if (this.format.getLabels() == true)
      {
        String str = getLabel(i);
        paramWriter.write(" '" + str + "'" + " \t");
      }
      float f;
      if (this.format.getWeights() == true)
      {
        f = getWeight(i);
        paramWriter.write(" " + f + " \t");
      }
      if (this.format.getConfidences() == true)
      {
        f = getConfidence(i);
        paramWriter.write(" " + f + " \t");
      }
      if (this.format.getIntervals() == true)
      {
        Interval localInterval = getInterval(i);
        if (localInterval == null) {
          paramWriter.write(" ()\t");
        } else {
          paramWriter.write(" " + localInterval.print() + "\t");
        }
      }
      paramWriter.write(" " + get(i) + ",\n");
    }
    paramWriter.write(";\n");
    paramWriter.write("END; [Splits]\n");
  }
  
//  public void read(NexusStreamParser paramNexusStreamParser, Taxa paramTaxa)
//    throws IOException
//  {
//    if (paramTaxa.getMustDetectLabels() == true) {
//      throw new IOException("line " + paramNexusStreamParser.lineno() + ": Can't read SPLITS block because no taxlabels given in TAXA block");
//    }
//    int[] arrayOfInt = new int[getNtax() + 1];
//    setNtax(paramTaxa.getNtax());
//    if (paramNexusStreamParser.peekMatchBeginBlock("st_splits")) {
//      paramNexusStreamParser.matchBeginBlock("st_splits");
//    } else {
//      paramNexusStreamParser.matchBeginBlock("Splits");
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("DIMENSIONS"))
//    {
//      paramNexusStreamParser.matchIgnoreCase("DIMENSIONS");
//      if (paramNexusStreamParser.peekMatchIgnoreCase("ntax=")) {
//        paramNexusStreamParser.matchIgnoreCase("ntax=" + getNtax());
//      }
//      if (paramNexusStreamParser.peekMatchIgnoreCase("nsplits="))
//      {
//        paramNexusStreamParser.matchIgnoreCase("nsplits=");
//        setNsplits(paramNexusStreamParser.getInt());
//      }
//      paramNexusStreamParser.matchIgnoreCase(";");
//    }
//    List localList;
//    if (paramNexusStreamParser.peekMatchIgnoreCase("FORMAT"))
//    {
//      localList = paramNexusStreamParser.getTokensLowerCase("format", ";");
//      this.format.labels = paramNexusStreamParser.findIgnoreCase(localList, "labels=no", false, this.format.labels);
//      this.format.labels = paramNexusStreamParser.findIgnoreCase(localList, "labels=left", true, this.format.labels);
//      this.format.weights = paramNexusStreamParser.findIgnoreCase(localList, "weights=no", false, this.format.weights);
//      this.format.weights = paramNexusStreamParser.findIgnoreCase(localList, "weights=yes", true, this.format.weights);
//      this.format.confidences = paramNexusStreamParser.findIgnoreCase(localList, "confidences=no", false, this.format.confidences);
//      this.format.confidences = paramNexusStreamParser.findIgnoreCase(localList, "confidences=yes", true, this.format.confidences);
//      this.format.intervals = paramNexusStreamParser.findIgnoreCase(localList, "intervals=no", false, this.format.intervals);
//      this.format.intervals = paramNexusStreamParser.findIgnoreCase(localList, "intervals=yes", true, this.format.intervals);
//      this.format.labels = paramNexusStreamParser.findIgnoreCase(localList, "no labels", false, this.format.labels);
//      this.format.labels = paramNexusStreamParser.findIgnoreCase(localList, "labels", true, this.format.labels);
//      this.format.weights = paramNexusStreamParser.findIgnoreCase(localList, "no weights", false, this.format.weights);
//      this.format.weights = paramNexusStreamParser.findIgnoreCase(localList, "weights", true, this.format.weights);
//      this.format.confidences = paramNexusStreamParser.findIgnoreCase(localList, "no confidences", false, this.format.confidences);
//      this.format.confidences = paramNexusStreamParser.findIgnoreCase(localList, "confidences", true, this.format.confidences);
//      this.format.intervals = paramNexusStreamParser.findIgnoreCase(localList, "no intervals", false, this.format.intervals);
//      this.format.intervals = paramNexusStreamParser.findIgnoreCase(localList, "intervals", true, this.format.intervals);
//      if (localList.size() != 0) {
//        throw new IOException("line " + paramNexusStreamParser.lineno() + ": `" + localList + "' unexpected in FORMAT");
//      }
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("threshold="))
//    {
//      paramNexusStreamParser.matchIgnoreCase("threshold=");
//      setThreshold((float)paramNexusStreamParser.getDouble());
//      paramNexusStreamParser.matchIgnoreCase(";");
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("PROPERTIES"))
//    {
//      localList = paramNexusStreamParser.getTokensLowerCase("properties", ";");
//      this.properties.fit = paramNexusStreamParser.findIgnoreCase(localList, "fit=", -1.0D, 100.0D, this.properties.fit);
//      if (paramNexusStreamParser.findIgnoreCase(localList, "weakly compatible", true, this.properties.compatibility == 3)) {
//        this.properties.compatibility = 3;
//      }
//      if (paramNexusStreamParser.findIgnoreCase(localList, "non compatible", true, this.properties.compatibility == 4)) {
//        this.properties.compatibility = 4;
//      }
//      if (paramNexusStreamParser.findIgnoreCase(localList, "compatible", true, this.properties.compatibility == 1)) {
//        this.properties.compatibility = 1;
//      }
//      if (paramNexusStreamParser.findIgnoreCase(localList, "cyclic", true, this.properties.compatibility == 2)) {
//        this.properties.compatibility = 2;
//      }
//      if (paramNexusStreamParser.findIgnoreCase(localList, "incompatible", true, this.properties.compatibility == 4)) {
//        this.properties.compatibility = 4;
//      }
//      if (paramNexusStreamParser.findIgnoreCase(localList, "leastsquares", true, getProperties().isLeastSquares())) {
//        getProperties().setLeastSquares(true);
//      } else {
//        getProperties().setLeastSquares(false);
//      }
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("CYCLE"))
//    {
//      paramNexusStreamParser.matchIgnoreCase("cycle");
//      arrayOfInt = new int[getNtax() + 1];
//      for (int i = 1; i <= getNtax(); i++) {
//        arrayOfInt[i] = paramNexusStreamParser.getInt();
//      }
//      paramNexusStreamParser.matchIgnoreCase(";");
//      try
//      {
//        setCycle(arrayOfInt);
//      }
//      catch (SplitsException localSplitsException)
//      {
//        //Basic.caught(localSplitsException);
//        //new Alert("Read-cycle failed: multiple occurrences of taxon");
//      }
//    }
//    if (paramNexusStreamParser.peekMatchIgnoreCase("matrix"))
//    {
//      paramNexusStreamParser.matchIgnoreCase("matrix");
//      readMatrix(paramNexusStreamParser);
//      paramNexusStreamParser.matchIgnoreCase(";");
//    }
//    paramNexusStreamParser.matchEndBlock();
//  }
//  
//  private void readMatrix(NexusStreamParser paramNexusStreamParser)
//    throws IOException
//  {
//    this.splits.setLabels(new Vector());
//    this.splits.setSplits(new Vector());
//    this.splits.setWeights(new Vector());
//    this.splits.setIntervals(new Vector());
//    int i = getNsplits();
//    setNsplits(0);
//    float f1 = 1.0F;
//    float f2 = -1.0F;
//    Interval localInterval = null;
//    for (int j = 1; j <= i; j++)
//    {
//      String str = null;
//      if (this.format.labels == true)
//      {
//        str = paramNexusStreamParser.getWordRespectCase();
//        if (str.equals("null")) {
//          str = null;
//        }
//      }
//      if (this.format.weights == true) {
//        f1 = (float)Math.max(0.0D, paramNexusStreamParser.getDouble());
//      }
//      if (this.format.confidences == true) {
//        f2 = (float)Math.max(0.0D, paramNexusStreamParser.getDouble());
//      }
//      if ((this.format.intervals == true) && (paramNexusStreamParser.peekMatchIgnoreCase("(")))
//      {
//        localInterval = new Interval();
//        paramNexusStreamParser.matchIgnoreCase("(");
//        localInterval.low = ((float)Math.max(0.0D, paramNexusStreamParser.getDouble()));
//        paramNexusStreamParser.matchIgnoreCase(",");
//        localInterval.high = ((float)Math.max(0.0D, paramNexusStreamParser.getDouble()));
//        paramNexusStreamParser.matchIgnoreCase(")");
//      }
//      TaxaSet localTaxaSet = new TaxaSet();
//      while (!paramNexusStreamParser.peekMatchIgnoreCase(","))
//      {
//        Integer localInteger = new Integer(paramNexusStreamParser.getWordRespectCase());
//        localTaxaSet.set(localInteger.intValue());
//      }
//      paramNexusStreamParser.matchIgnoreCase(",");
//      if ((localTaxaSet.cardinality() == 0) || (localTaxaSet.cardinality() == getNtax())) {
//        throw new IOException("line " + paramNexusStreamParser.lineno() + ": non-split of size " + localTaxaSet.cardinality());
//      }
//      if (f2 == -1.0F) {
//        f2 = 1.0F;
//      }
//      add(localTaxaSet, f1, f2, localInterval, str);
//    }
//  }
//  
  public String toString()
  {
    return new String("[Splits, ntax=" + getNtax() + ", nsplits=" + getNsplits() + "]");
  }
  
  public String toString(Taxa paramTaxa)
  {
    StringWriter localStringWriter = new StringWriter();
    try
    {
      write(localStringWriter, paramTaxa);
    }
    catch (IOException localIOException)
    {
      //Basic.caught(localIOException);
    }
    return localStringWriter.toString();
  }
  
//  public Splits clone(Taxa paramTaxa)
//  {
//    Splits localSplits = new Splits();
//    if ((getNtax() > 0) && (getNsplits() > 0)) {
//      localSplits.copy(paramTaxa, this);
//    }
//    return localSplits;
//  }
  
//  public void copy(Taxa paramTaxa, Splits paramSplits)
//  {
//    boolean bool1 = paramSplits.getFormat().getConfidences();
//    boolean bool2 = paramSplits.getFormat().getWeights();
//    boolean bool3 = paramSplits.getFormat().getLabels();
//    paramSplits.getFormat().setWeights(true);
//    paramSplits.getFormat().setConfidences(true);
//    paramSplits.getFormat().setLabels(true);
//    try
//    {
//      StringWriter localStringWriter = new StringWriter();
//      paramSplits.write(localStringWriter, paramTaxa);
//      read(new NexusStreamParser(new StringReader(localStringWriter.toString())), paramTaxa);
//    }
//    catch (Exception localException)
//    {
//      //Basic.caught(localException);
//    }
//    getFormat().setConfidences(bool1);
//    getFormat().setWeights(bool2);
//    getFormat().setLabels(bool3);
//    paramSplits.getFormat().setConfidences(bool1);
//    paramSplits.getFormat().setWeights(bool2);
//    paramSplits.getFormat().setLabels(bool3);
//  }
//  
  public boolean getFormatSwitchValue(String paramString)
  {
    if (paramString.equalsIgnoreCase("labels")) {
      return getFormat().getLabels();
    }
    if (paramString.equalsIgnoreCase("weights")) {
      return getFormat().getWeights();
    }
    if (paramString.equalsIgnoreCase("confidences")) {
      return getFormat().getConfidences();
    }
    return true;
  }
  
//  public void hideTaxa(Taxa paramTaxa, TaxaSet paramTaxaSet)
//  {
//    if (((paramTaxaSet == null) || (paramTaxaSet.cardinality() == 0)) && (this.originalSplits == null)) {
//      return;
//    }
//    if (this.originalSplits == null) {
//      this.originalSplits = clone(paramTaxa);
//    }
//    Taxa localTaxa = Taxa.getInduced(paramTaxa, paramTaxaSet);
//    if ((this.previousTaxa != null) && (localTaxa.equals(this.previousTaxa))) {
//      return;
//    }
//    this.previousTaxa = localTaxa;
//    int[] arrayOfInt1 = new int[paramTaxa.getNtax() + 1];
//    int i = 0;
//    for (int j = 1; j <= paramTaxa.getNtax(); j++) {
//      if ((paramTaxaSet == null) || (!paramTaxaSet.get(j))) {
//        arrayOfInt1[j] = (++i);
//      }
//    }
//    int[] arrayOfInt2 = null;
//    int[] arrayOfInt3 = getCycle();
//    if (getCycle() != null)
//    {
//      arrayOfInt2 = new int[localTaxa.getNtax() + 1];
//      int k = 0;
//      for (int m = 1; m <= getNtax(); m++)
//      {
//        int n = arrayOfInt3[m];
//        if (arrayOfInt1[n] != 0) {
//          arrayOfInt2[(++k)] = arrayOfInt1[n];
//        }
//      }
//    }
//    clear();
//    setNtax(localTaxa.getNtax());
//    HashMap localHashMap = new HashMap();
//    int[] arrayOfInt4 = new int[this.originalSplits.getNsplits() + 1];
//    int n = 0;
//    for (int i1 = 1; i1 <= this.originalSplits.getNsplits(); i1++)
//    {
//      TaxaSet localTaxaSet = new TaxaSet();
//      for (int i3 = 1; i3 <= paramTaxa.getNtax(); i3++) {
//        if (((paramTaxaSet == null) || (!paramTaxaSet.get(i3))) && (this.originalSplits.get(i1).get(i3)) && (arrayOfInt1[i3] != 0)) {
//          localTaxaSet.set(arrayOfInt1[i3]);
//        }
//      }
//      if (!localTaxaSet.get(1)) {
//        localTaxaSet = localTaxaSet.getComplement(localTaxa.getNtax());
//      }
//      if ((localTaxaSet.cardinality() > 0) && (localTaxaSet.cardinality() < localTaxa.getNtax())) {
//        if (localHashMap.containsKey(localTaxaSet))
//        {
//          Integer i3 = ((Integer)localHashMap.get(localTaxaSet)).intValue();
//          setWeight(i3, getWeight(i3) + this.originalSplits.getWeight(i1));
//          setConfidence(i3, (getConfidence(i3) + this.originalSplits.getConfidence(i1)) / 2.0F);
//          if ((getLabel(i3) != null) && (this.originalSplits.getLabel(i1) != null)) {
//            setLabel(i3, getLabel(i3) + "_" + this.originalSplits.getLabel(i1));
//          }
//          arrayOfInt4[i3] += 1;
//        }
//        else
//        {
//          add(localTaxaSet, this.originalSplits.getWeight(i1), this.originalSplits.getConfidence(i1), this.originalSplits.getLabel(i1));
//          localHashMap.put(localTaxaSet, new Integer(++n));
//          arrayOfInt4[n] = 1;
//        }
//      }
//    }
//    try
//    {
//      setCycle(arrayOfInt2);
//    }
//    catch (SplitsException localSplitsException)
//    {
//      //Basic.caught(localSplitsException);
//    }
//    getProperties().setFit(-1.0D);
//    getProperties().setLSFit(-1.0D);
//    for (int i2 = 1; i2 <= getNsplits(); i2++) {
//      setConfidence(i2, getConfidence(i2) / arrayOfInt4[i2]);
//    }
//  }
//  
//  public void restoreOriginal(Taxa paramTaxa)
//  {
//    copy(paramTaxa, this.originalSplits);
//    this.previousTaxa = paramTaxa;
//  }
//  
//  public void setOriginal(Taxa paramTaxa)
//  {
//    this.originalSplits = clone(paramTaxa);
//    this.previousTaxa = null;
//  }
//  
  public Splits getOriginal()
  {
    return this.originalSplits;
  }
  
//  public void hideSplits(Taxa paramTaxa, BitSet paramBitSet)
//  {
//    if (paramBitSet.cardinality() > 0)
//    {
//      if (this.originalSplits == null) {
//        setOriginal(paramTaxa);
//      }
//      Splits localSplits = clone(paramTaxa);
//      clear();
//      setNtax(paramTaxa.getNtax());
//      for (int i = 1; i <= localSplits.getNsplits(); i++) {
//        if (!paramBitSet.get(i)) {
//          add(localSplits.get(i), localSplits.getWeight(i), localSplits.getConfidence(i), localSplits.getLabel(i));
//        }
//      }
//    }
//  }
  
  public class Properties
    implements Cloneable
  {
    public static final int COMPATIBLE = 1;
    public static final int CYCLIC = 2;
    public static final int WEAKLY_COMPATIBLE = 3;
    public static final int INCOMPATIBLE = 4;
    public static final int UNKNOWN = 5;
    int compatibility = 5;
    private double fit = -1.0D;
    private double lsfit = -1.0D;
    private boolean leastSquares = false;
    
    public boolean isLeastSquares()
    {
      return this.leastSquares;
    }
    
    public void setLeastSquares(boolean paramBoolean)
    {
      this.leastSquares = paramBoolean;
    }
    
    public Properties() {}
    
    public double getFit()
    {
      return this.fit;
    }
    
    public void setFit(double paramDouble)
    {
      this.fit = paramDouble;
    }
    
    public double getLSFit()
    {
      return this.lsfit;
    }
    
    public void setLSFit(double paramDouble)
    {
      this.lsfit = paramDouble;
    }
    
    public int getCompatibility()
    {
      return this.compatibility;
    }
    
    public void setCompatibility(int paramInt)
    {
      this.compatibility = paramInt;
    }
  }
  
  public class Format
    implements Cloneable
  {
    private boolean labels = false;
    private boolean weights = true;
    private boolean confidences = false;
    private boolean intervals = false;
    
    public Format() {}
    
    public boolean getLabels()
    {
      return this.labels;
    }
    
    public boolean getWeights()
    {
      return this.weights;
    }
    
    public void setLabels(boolean paramBoolean)
    {
      this.labels = paramBoolean;
    }
    
    public void setWeights(boolean paramBoolean)
    {
      this.weights = paramBoolean;
    }
    
    public boolean getConfidences()
    {
      return this.confidences;
    }
    
    public void setConfidences(boolean paramBoolean)
    {
      this.confidences = paramBoolean;
    }
    
    public boolean getIntervals()
    {
      return this.intervals;
    }
    
    public void setIntervals(boolean paramBoolean)
    {
      this.intervals = paramBoolean;
    }
  }
}
