# Part 2: Data Files

---

## Getting Ready
```
cd ~/Developer/KorfLab/
git clone https://github.com/iankorf/E.coli
```

Since there's only one file for each type of file, wildcards were used instead of full file names.

## Question 1

Not knowing if the file has duplicates, my initial code:

```
gunzip -c *.faa.gz | grep "^>" | sort | uniq | wc -l
```

Output:

```
4298
```

Doublechecking:

```
gunzip -c *.faa.gz | grep -c "^>"
```

Output:

```
4298
```

So the file has no duplicates!

**Answer:**

According to the `*.faa` file, it encodes 4298 proteins.

## Question 2

My code:

```
gunzip -c *.gbff.gz | grep -c "^     CDS"
```

Output:

```
4315
```

**Answer:**

According to the `*.gbff` file, there are 4315 coding sequences?

## Question 3

**Answer:**

There are 4315 coding sequences but it encodes only 4298 proteins.

This difference in the two may have been the result of pseudogene, which are gene-looking sequences that would not contribute to the protein count, hence a smaller number of encoded proteins than coding sequences. 

Also, alternative splicing might have led different coding sequences to the same protein. 

Moreover, though rare, there's always a chance for annotation errors that might contributed to this discrepancy.

## Question 4

My code:

```
gunzip -c *.gff.gz | grep -c "tRNA"
```

Output:

```
332
```

**Answer:**

There are 332 tRNA features found in the `*.gff` file.

## Question 5

My code:

```
gunzip -c *.gff.gz | grep -v "^#" | cut -f3 | sort | uniq -c
```

Output:

```
4337 CDS
 207 exon
4494 gene
  50 mobile_genetic_element
  99 ncRNA
   1 origin_of_replication
 145 pseudogene
  22 rRNA
   1 region
  48 sequence_feature
  86 tRNA
```

**Answer:**

| Feature Type           | Count |
| ---------------------- | ----- |
| CDS                    | 4337  |
| exon                   | 207   |
| gene                   | 4494  |
| mobile_genetic_element | 50    |
| ncRNA                  | 99    |
| origin_of_replication  | 1     |
| pseudogene             | 145   |
| rRNA                   | 22    |
| region                 | 1     |
| sequence_feature       | 48    |
| tRNA                   | 86    |

---