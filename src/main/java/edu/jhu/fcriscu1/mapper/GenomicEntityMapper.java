package org.nygenome.genomic.tool.mapper;

import com.github.davidmoten.rtree.Entry;
import com.github.davidmoten.rtree.RTree;
import com.github.davidmoten.rtree.geometry.Geometries;
import com.github.davidmoten.rtree.geometry.Rectangle;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.google.common.collect.Maps;
import io.reactivex.Observable;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import lombok.extern.log4j.Log4j;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.validator.routines.DoubleValidator;
import org.apache.commons.validator.routines.FloatValidator;
import org.eclipse.collections.impl.factory.Lists;
import scala.Tuple2;

@Log4j
public class GenomicEntityMapper implements IdMapper {

  private Map<String, String> entrezMap = Suppliers.memoize(new EntrezIDSupplier()).get();
  private Map<String, Tuple2<String, String>>
      gnMap = Suppliers.memoize(new EnsemblNameMapSupplier()).get();
  private RTree<String, Rectangle> genePositionTree = Suppliers.memoize(new GenomicNameSupplier())
      .get();

  public static void main(String... args) {
    GenomicEntityMapper test = new GenomicEntityMapper();

    log.info(test.findGeneNameByGenomicPosition("1", "42817156", "1")); //ERMAP
    log.info(test.findGeneNameByGenomicPosition("1", "42817156", "-1")); //CCDC23
    log.info(test.findGeneNameByGenomicPosition("8", "13566900", "1")); //C8orf48
    log.info(test.findGeneNameByGenomicPosition("8", "1000", "1")); // intergenic
    log.info(test.findGeneNameByGenomicPosition("X", "91982660", "-1")); //VDAC1P3
    log.info(test.findGeneNameByGenomicPosition("Y", "6246223", "1")); //TSPY2
    log.info(test.findGeneNameByGenomicPosition("17", "43044296", "-1")); //BRCA1
    // test for genomic range
    //43044295	43170245
    Stopwatch stopwatch = Stopwatch.createStarted();
    test.findGeneByGenomicRange("17","43044295","43170245")
        .forEach(gene -> log.info("Gene: " +gene));
    log.info("Genomic range search required: " +stopwatch.elapsed(TimeUnit.MILLISECONDS) +" milliseconds");
  }


  private Function<String, Double> resolveChromosomeNumberFunction = (chromosome) -> {
    if (chromosome.toUpperCase().equals("X")) {
      return 23.0D;
    }
    if (chromosome.toUpperCase().equals("Y")) {
      return 24.0D;
    }
    return DoubleValidator.getInstance().validate(chromosome);
  };

  public List<String> findGeneByGenomicRange(String chromosome, String startPosition, String endPosition){
    // TODO: validate parameters
    List<String> geneList = Lists.mutable.empty();
    Double x = resolveChromosomeNumberFunction.apply(chromosome);
    Double y1 = DoubleValidator.getInstance().validate(startPosition);
    Double y2 = DoubleValidator.getInstance().validate(endPosition);
    this.genePositionTree.search(Geometries.rectangle(x,y1,x,y2))
        .doOnNext(entry-> geneList.add(entry.value()))
        .doOnError(e -> log.error(e.getMessage()))
        .subscribe();
    // look for genes on minys strand

    this.genePositionTree.search(Geometries.rectangle(x,(y2*-1.0D),x,(y1*-1.0D)))
        .doOnNext(entry-> geneList.add(entry.value()))
        .doOnError(e -> log.error(e.getMessage()))
        .subscribe();
   log.info("Search by genomic range found " + geneList.size() +"  genes on pkus and minus strands");
    return  geneList;
  }

  public String findGeneNameByGenomicPosition(String chromosome, String position, String strand) {
    Preconditions
        .checkArgument(!Strings.isNullOrEmpty(chromosome), "A chromosome name is required");
    Preconditions
        .checkArgument(!Strings.isNullOrEmpty(position), "A chromosome position is required");
    if (Strings.isNullOrEmpty(strand)) {
      strand = "1";
    }
    Double x = resolveChromosomeNumberFunction.apply(chromosome);
    // negate the position if it relates to the negative strand
    Double y = (strand.equals("1")) ? DoubleValidator.getInstance().validate(position) :
        DoubleValidator.getInstance().validate(position) * -1.0D;
    if (null == x || null == y) {
      log.error("Invalid chromosome or position was provided");
      return "";
    }

    Entry<String, Rectangle> result = null;
    try {
      result = this.genePositionTree.search(Geometries.point(x, y))
          .toBlocking().first();
    } catch (Exception e) {
      return CommonNames.INTERGENIC; // a gene name wasn't found
    }
    return result.value();
  }

  @Override
  public String symbolToEntrezID(String geneSymbol) {
    if (!Strings.isNullOrEmpty(geneSymbol)) {
      return this.entrezMap.getOrDefault(geneSymbol, "");
    }
    return "";
  }

  @Override
  public String entrezIDToSymbol(String entrezID) throws Exception {
    Preconditions.checkArgument(!Strings.isNullOrEmpty(entrezID));
    for (Map.Entry<String, String> entry : this.entrezMap.entrySet()) {
      if (entry.getValue().equals(entrezID)) {
        return entry.getKey();
      }
    }
    return "";
  }

    /*
    method to return the HUGO Symbol and EntrezID for a specified Ensembl ID
    return object is a Tuple2 containing the HUGO Gene Symbol & the EntrezID
    an empty tuple is returned if a mapping cannot be completed
     */

  public Tuple2<String, String> ensemblToHugoSymbolAndEntrezID(String ensemblID) {
    Preconditions.checkArgument(!Strings.isNullOrEmpty(ensemblID));
    return (this.gnMap.containsKey(ensemblID)) ? this.gnMap.get(ensemblID) :
        new Tuple2<String, String>("", "");
  }

  /*
inner class responsible for creating an RTree based on a file of genomic entities and
their chromosome name, position, and strand
 */
  public class GenomicNameSupplier implements Supplier<RTree<String, Rectangle>> {

    private InputStreamReader reader;
    private final String ENSEMBL_GENE_FILE = "/ensembl_grch38_gene.tsv";
    private RTree<String, Rectangle> genePositionTree;

    public GenomicNameSupplier() {
      reader = new InputStreamReader((this.getClass().getResourceAsStream(ENSEMBL_GENE_FILE)));
      this.genePositionTree = RTree.create();
    }

    private Predicate<String> validChromosomePredicate =
        CommonNames.validChromosomeSet::contains;

    private Predicate <String> validStringPredicate = (s) ->
      !Strings.isNullOrEmpty(s) && s.length()>1;

//Ensembl Gene ID	Chromosome	Gene start	Gene end	HGNC symbol	Strand
    private Function<CSVRecord, Tuple2<String, Optional<Rectangle>>> resolveRectangleFunction
        = (record) -> new Tuple2<>(record.get("HGNC symbol"),
        resolveGenePosition(record.get("Chromosome"),
            record.get("Gene start"), record.get("Gene end"), record.get("Strand"))
    );

    /*
    Private Consumer to add a Rectangle to the Tree if present
     */
    private Consumer<Tuple2<String, Optional<Rectangle>>> rectangleTupleConsumer = nameRectTuple ->
        nameRectTuple._2().ifPresent(rect ->
            genePositionTree = genePositionTree.add(nameRectTuple._1(), rect));

    @Override
    public RTree<String, Rectangle> get() {
      try {
        final CSVParser parser = new CSVParser(this.reader, CSVFormat.TDF.withHeader());
        Observable.fromIterable(parser.getRecords())
            .filter(record -> validChromosomePredicate.test(record.get("Chromosome")))
            .filter(record -> validStringPredicate.test(record.get("HGNC symbol")))
            // transform chromosome coordinates to RTree Rectangles
            .map(record -> resolveRectangleFunction.apply(record))
            .doOnNext(rectangleTupleConsumer::accept)
            .doOnError(throwable -> log.error(throwable.getMessage()))
            .doOnComplete(() -> log.info("Completed RTree, size = " + genePositionTree.size()))
            .subscribe();
      } catch (IOException e) {
        log.error(e.getMessage());
        e.printStackTrace();
      }
      return this.genePositionTree;
    }

    /*
    private method to create a rectangle with zero height representing the "area" covered by
    the genomic entity - really just the length
    n.b. starting and end coordinates are reversed for features on the minus strand
     */
    private Optional<Rectangle> resolveGenePosition(String chromosome, String startPos,
        String endPos, String strand) {
      if (Strings.isNullOrEmpty(strand)) {
        strand = "1";
      }
      Float s = FloatValidator.getInstance().validate(strand);
      Float y1;
      Float y2;
      if (s > 0.0) {
        y1 = FloatValidator.getInstance().validate(startPos);
        y2 = FloatValidator.getInstance().validate(endPos);
      } else {
        y1 = FloatValidator.getInstance().validate(endPos) * s;
        y2 = FloatValidator.getInstance().validate(startPos) * s;
      }

      // convert  X & Y chromosome to numeric values (X=22, Y= 23)
      if (chromosome.toUpperCase().equals("X")) {
        chromosome = "23";
      } else if (chromosome.toUpperCase().equals("Y")) {
        chromosome = "24";
      }
      Float x1 = FloatValidator.getInstance().validate(chromosome);

      if (null == x1 || null == y1 || null == y2) {
        return Optional.empty();
      }
      Float x2 = x1; // pseudo rectangle
      return Optional.of(Geometries.rectangle(x1, y1, x2, y2));
    }
  }

  public class EntrezIDSupplier implements Supplier<Map<String, String>> {

    private InputStreamReader reader;
    private final String ESEMBL_GENE_LIST_FILE = "ensembl_grch38_gene.tsv";

    public EntrezIDSupplier() {
      ClassLoader classLoader = getClass().getClassLoader();
      //File file = new File(classLoader.getResource(HUGO_GENE_FILE).getFile());
      reader = new InputStreamReader(
          GenomicEntityMapper.class.getClassLoader().getResourceAsStream(ESEMBL_GENE_LIST_FILE));
    }

    /*
     public method to create and supply a Map of HUGO Symbols keys
     and  Entrez id as values
     n.b. the Entrez ID is treated as a numeric String
     */
    @Override
    public Map<String, String> get() {

      Map<String, String> hugoMap = Maps.newHashMap();
      try {
        final CSVParser parser = new CSVParser(this.reader, CSVFormat.TDF.withHeader());
        for (CSVRecord record : parser) {

          hugoMap.put(record.get(0), record.get(1));
        }

      } catch (IOException e) {
        log.error(e.getMessage());
        e.printStackTrace();
      }
      return hugoMap;
    }
  }

  public class EnsemblNameMapSupplier implements Supplier<Map<String, Tuple2<String, String>>> {

    private static final String ENSEMBL_GENE_FILE = "/HGNC_Ensembl.tsv";
    private InputStreamReader reader;

    public EnsemblNameMapSupplier() {
      reader = new InputStreamReader((this.getClass().getResourceAsStream(ENSEMBL_GENE_FILE)));
    }

    /*
    public method to create and supply a Map of Ensemble gene ids as keys
    and Tuple2s containing Hugo symbol & Entrez id as values
    */
    @Override
    public Map<String, Tuple2<String, String>> get() {

      Map<String, Tuple2<String, String>> ensemblMap = Maps.newHashMap();
      try {
        final CSVParser parser = new CSVParser(this.reader, CSVFormat.TDF.withHeader());
        for (CSVRecord record : parser) {
          ensemblMap
              .put(record.get("Ensembl"), new Tuple2(record.get("Symbol"), record.get("Entrez")));
        }

      } catch (IOException e) {
        log.error(e.getMessage());
        e.printStackTrace();
      }

      return ensemblMap;

    }
  }



}
