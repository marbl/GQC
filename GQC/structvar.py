import sys
import logging
import pybedtools
from GQC import bedtoolslib

logger = logging.getLogger(__name__)

# reminders: in aligndata, (1) all coordinates are 1-based, (2) strand is "+" or "-", (3) querystart is the query's lower coordinate,
# so doesn't correspond to targetstart if alignment is on the reverse strand

def write_structural_errors(aligndata:list, refobj, queryobj, outputdict, bedobjectdict, bmstats, args)->str:

    aligndict = {}
    current_align = None
    structvarbedstring = ""
    for align in sorted(aligndata, key=lambda a: (a["target"], a["targetstart"], a["targetend"])):
        refentry = align["target"]
        if refentry not in aligndict:
            aligndict[refentry] = [align]
        else:
            aligndict[refentry].append(align)
        query = align["query"]
        refstart = align["targetstart"]
        refend = align["targetend"]
        querystart = align["querystart"]
        queryend = align["queryend"]
        strand = align["strand"]
        if current_align is not None:
            ref1 = current_align["targetend"]
            ref2 = refstart
            refdiff = ref2 - ref1
            if refentry == current_align["target"] and query == current_align["query"] and strand == current_align["strand"]:
                if strand == "+":
                    query1 = current_align["queryend"]
                    query2 = querystart
                    querydiff = query2 - query1
                    netdiff = querydiff - refdiff
                else:
                    querydiff = queryend - current_align["querystart"]
                    netdiff = -1.0 * querydiff - refdiff
                    query1 = querystart
                    query2 = current_align["queryend"]
 
                if query1 <= query2:
                    minquery = query1
                    maxquery = query2
                else:
                    minquery = query2
                    maxquery = query1

                # check for overlap with N's in the test assembly, and if there are none, add to "same contig" structural variants:
                queryspanobj = pybedtools.BedTool(query + "\t" + str(minquery) + "\t" + str(maxquery) + "\n", from_string = True)
                overlaptestnobj = bedtoolslib.intersectintervals(queryspanobj, bedobjectdict["testnregions"])

                if overlaptestnobj.count == 0:
                    if refdiff < querydiff: # refdiff less than querydiff (insertion), netshift positive
                        if refdiff > 0: # print ref1 before ref2
                            structvarbedstring = structvarbedstring + refentry + "\t" + str(ref1 - 1) + "\t" + str(ref2) + "\tSameContigInsertion\t" + query + "\t" + str(minquery) + "\t" + str(maxquery) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n"
                        else:
                            structvarbedstring = structvarbedstring + refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tSameContigInsertion\t" + query + "\t" + str(minquery) + "\t" + str(maxquery) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n"
                    else: # refdiff greater than than querydiff (deletion), netshift negative
                        if refdiff > 0:
                            structvarbedstring = structvarbedstring + refentry + "\t" + str(current_align["targetend"] - 1) + "\t" + str(refstart) + "\tSameContigDeletion\t" + query + "\t" + str(minquery) + "\t" + str(maxquery) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n"
                        else:
                            structvarbedstring = structvarbedstring + refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tSameContigDeletion\t" + query + "\t" + str(minquery) + "\t" + str(maxquery) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n"
            elif refentry == current_align["target"]: # strand switch or new contig:
                queryentries = query + "/" + current_align["query"]
                strands = strand + "/" + current_align["strand"]
                if refdiff > 0:
                    structvarbedstring = structvarbedstring + refentry + "\t" + str(current_align["targetend"] - 1) + "\t" + str(refstart) + "\tBetweenContigDeletion\t" + queryentries + "\t.\t.\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\tNA\t" + strands + "\n"
                else:
                    structvarbedstring = structvarbedstring + refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tBetweenContigInsertion\t" + queryentries + "\t.\t.\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\tNA\t" +strands + "\n"
        current_align = align

    allstructvarbedobject = pybedtools.BedTool(structvarbedstring, from_string = True).sort()
    # don't report structural variants that overlap excluded regions or assembly N's
    nonexcludedstructvarbedobject = bedtoolslib.intersectintervals(allstructvarbedobject, bedobjectdict["allexcludedregions"], v=True, wa=True)
    nonexcludedstructvarbedobject.saveas(outputdict["structvariantbed"])
    # Can't intersect with test assembly coordinates--this is in benchmark space!
    #nonexcludednonnstructvarbedobject = bedtoolslib.intersectintervals(nonexcludedstructvarbedobject, bedobjectdict["testnregions"], v=True, wa=True)
    #nonexcludednonnstructvarbedobject.saveas(outputdict["structvariantbed"])

    return 0

