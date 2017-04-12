/*
 *  Copyright (c) 2014 Memorial Sloan-Kettering Cancer Center.
 * 
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 *  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 *  documentation provided hereunder is on an "as is" basis, and
 *  Memorial Sloan-Kettering Cancer Center 
 *  has no obligations to provide maintenance, support,
 *  updates, enhancements or modifications.  In no event shall
 *  Memorial Sloan-Kettering Cancer Center
 *  be liable to any party for direct, indirect, special,
 *  incidental or consequential damages, including lost profits, arising
 *  out of the use of this software and its documentation, even if
 *  Memorial Sloan-Kettering Cancer Center 
 *  has been advised of the possibility of such damage.
 */
package edu.jhu.fcriscu1.util;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.*;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * represents a collection of common values
 * @author fcriscuolo
 */
public interface StagingCommonNames {

    public static final String HUGO_COLUMNNAME = "Hugo_symbol";
    public static final String INTERGENIC = "intergenic";

    // common file extensions
    public static final String TEXT_FILE_EXTENSION = "txt";
    public static final String EXCEL_FILE_EXTENSION = "xls";
    public static final String EXCELX_FILE_EXTENSION = "xlsx";
    public static final String TSV_FILE_EXTENSION = "tsv";
    public static final String XML_FILE_EXTENSION = "xml";
    public static final String COMPRESSED_FILE_EXTENSION = "gz";

    // support for legacy file extension use
    public static final String stagingFileExtension = TEXT_FILE_EXTENSION;
    public static final String xmlExtension = XML_FILE_EXTENSION;


    // COMMON SPLITTERS & JOINERS
    public static final Splitter tabSplitter = Splitter.on("\t");
    public static final Splitter lineSplitter = Splitter.on("\n").trimResults();
    public static final Splitter dashSplitter = Splitter.on("-").trimResults();
    public static final Splitter blankSplitter = Splitter.on(" ");
    public static final Joiner tabJoiner = Joiner.on('\t').useForNull(" ");
    public static final Joiner commaJoiner = Joiner.on(',').useForNull(" ");
    public static final Splitter commaSplitter= Splitter.on(',');
    public static final Joiner blankJoiner = Joiner.on(" ");
    public static final Joiner dashJoiner = Joiner.on("-");
    public static final Joiner lineJoiner = Joiner.on("\n");
    public static final Splitter posSplitter = Splitter.on(':');
    public static final Splitter semicolonSplitter = Splitter.on(';');
    public final Joiner pathJoiner =
            Joiner.on(System.getProperty("file.separator"));
    public final Splitter pathSplitter =
            Splitter.on(System.getProperty("file.separator"));
    public final Pattern tabPattern = Pattern.compile("\t");


    // length of human chromosomes
    // http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
    public static final ImmutableMap<String, Long> chromosomeLengthMap = new ImmutableMap.Builder<String, Long>()
            .put("1", Long.valueOf(248_956_422))
                    .put("2", Long.valueOf(242_193_529))
                    .put("3", Long.valueOf(198_295_559))
                    .put("4", Long.valueOf(190_214_555))
                    .put("5", Long.valueOf(181_538_259))
                    .put("6", Long.valueOf(170_805_979))
                    .put("7", Long.valueOf(159_345_973))
                    .put("8", Long.valueOf(145_138_636))
                    .put("9", Long.valueOf(138_394_717))
                    .put("10", Long.valueOf(133_797_42))
                    .put("11", Long.valueOf(135_086_622))  // this is correct
                    .put("12", Long.valueOf(133_75_309))
                    .put("13", Long.valueOf(114_364_328))
                    .put("14", Long.valueOf(107_043_718))
                    .put("15", Long.valueOf(101_991_189))
                    .put("16", Long.valueOf(90_338_345))
                    .put("17", Long.valueOf(83_257_441))
                    .put("18", Long.valueOf(80_373_285))
                    .put("19", Long.valueOf(58_617_616))
                    .put("20", Long.valueOf(64_444_167))
                    .put("21", Long.valueOf(46_709_983))
                    .put("22", Long.valueOf(50_818_468))
                    .put("X", Long.valueOf(156_040_895))  // upper and lower case for x & y
                    .put("Y", Long.valueOf(57_227_415))
                    .put("x", Long.valueOf(156_040_895))
                    .put("y", Long.valueOf(57_227_415)).build();
    // valid chromosome values
    public static final Set<String> validChromosomeSet = Sets.newHashSet("1","2","3","4","5","6",
           "7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y",
           "x","y" );
}