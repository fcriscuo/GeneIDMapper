package org.nygenome.genomic.tool.mapper;

import scala.Tuple2;

public interface IdMapper {
  /**
   *
   * @param chromosome
   * @param position
   * @param strand
   * @return
   */

  public String findGeneNameByGenomicPosition(String chromosome, String position, String strand);

  /**
   * For the given symbol, return id.
   *
   * @param geneSymbol String
   * @return String
   * @throws Exception
   */
  String symbolToEntrezID(String geneSymbol)  throws Exception;

  /**
   * For the entrezID, return symbol.
   *
   * @param entrezID String
   * @return String
   * @throws Exception
   */
  String entrezIDToSymbol(String entrezID) throws Exception;

  /**
   * returns the Gene Symbol and Entrez ID for a specified Ensembl ID
   * @param ensemblID
   * @return
   */

  Tuple2<String,String> ensemblToHugoSymbolAndEntrezID(String ensemblID);
}
