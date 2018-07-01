cat $1 | grep -v \> | grep -E -o "G|C" | wc -l > /tmp/GC_content
cat $1 | grep -v \> | grep -E -o "G|C|A|T" | wc -l >> /tmp/GC_content
cat /tmp/GC_content | awk '{getline total; print $0/total}'
