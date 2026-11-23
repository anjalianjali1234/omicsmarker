# Compare runs in MongoDB

In the omicsmarker database, aggregate folds within each run:

```js
db.runs.aggregate([
  { $unwind: "$metrics" },
  { $group: { _id: {run: "$_id", target: "$target", model: "$metrics.model", assay: "$metrics.features"}, auroc: {$avg: "$metrics.auroc"} } },
  { $sort: { auroc: -1 } }
])
```
